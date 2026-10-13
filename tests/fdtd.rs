#![cfg(not(target_arch = "wasm32"))]

use fieldgrid::Gpu;
use fieldgrid::fdtd::{
    self, C0, Debye, Edge, Element, Extras, Grid, Lumped, Materials, Media, PortDef, Pulse, Record,
    Resistive, Sheets, Sim, idx,
};
use fieldgrid::gpu::{Options, Software};
use std::sync::Mutex;

static DEVICE: Mutex<()> = Mutex::new(());

fn gpu() -> Option<Gpu> {
    let options = Options { software: Software::Allow, storage_buffers: 12, ..Options::default() };
    match Gpu::open_blocking(&options) {
        Ok(g) => Some(g),
        Err(e) if !cfg!(feature = "gpu-tests") => {
            eprintln!("skipped: {e}");
            None
        }
        Err(e) => panic!("{e}"),
    }
}

fn graded(n: usize, start: f64, grow: f64) -> Vec<f64> {
    let mut out = vec![0.0];
    let mut d = start;
    for _ in 1..n {
        out.push(out.last().unwrap() + d);
        d *= grow;
    }
    out
}

fn sheet_table() -> Sheets {
    let table = vec![
        1.0, 1.0, 2.0, 0.0, 0.9, 0.05, 0.8, 0.1, 0.02, 0.01, 1.0, 0.2, 0.1, 0.05, 0.02, 1.2, 0.1,
        0.08, 0.03, 0.01,
    ];
    Sheets { table, states: 3 }
}

fn everything() -> Sim {
    let grid = Grid {
        x: graded(26, 1e-3, 1.01),
        y: graded(22, 1e-3, 1.0),
        z: graded(20, 0.8e-3, 1.02),
        pml: 4,
    };
    let n = grid.dims();
    let mut mats = Materials::new(&grid);
    let cells = [n[0] - 1, n[1] - 1, n[2] - 1];
    for i in 0..cells[0] {
        for j in 0..cells[1] {
            for k in 6..9 {
                let id = (i * cells[1] + j) * cells[2] + k;
                mats.eps[id] = 4.0;
                mats.inf[id] = 3.6;
                mats.debye[id] = 0.4;
                mats.sigma[id] = 0.002;
            }
        }
    }
    let debye = Debye::band(1e9, 2e10);
    let pec = |c: usize, p: [usize; 3]| c != 2 && p[2] == 6 && p[0] < 7;
    let lumped = vec![
        Lumped {
            name: "c".into(),
            edges: vec![Edge { comp: 2, at: [14, 11, 6] }],
            element: Element::Capacitor(1e-12),
        },
        Lumped {
            name: "r".into(),
            edges: vec![Edge { comp: 2, at: [16, 11, 6] }],
            element: Element::Resistor(30.0),
        },
        Lumped {
            name: "l".into(),
            edges: vec![Edge { comp: 2, at: [18, 11, 6] }],
            element: Element::Inductor(2e-9),
        },
    ];
    let column =
        |i: usize, j: usize| (6..9).map(|k| Edge { comp: 2, at: [i, j, k] }).collect::<Vec<_>>();
    let ports = vec![
        PortDef {
            name: "a".into(),
            columns: vec![column(10, 10), column(10, 11)],
            r: 50.0,
            reference_above: false,
        },
        PortDef { name: "b".into(), columns: vec![column(20, 12)], r: 75.0, reference_above: true },
    ];
    let resistive: Vec<Resistive> = (8..14)
        .map(|i| Resistive {
            comp: 0,
            id: idx(n, i, 8, 9),
            r: 0.4,
            adaptive: i % 2 == 0,
            faces: [0, 1],
        })
        .collect();
    let media = Media { sheets: sheet_table(), debye };
    Sim::new(grid, &mats, &pec, &lumped, ports, &resistive, media)
}

fn extras(sim: &Sim, steps: usize, fused_e: bool) -> Extras {
    Extras {
        max_steps: steps,
        decay_db: f64::INFINITY,
        freqs: vec![3e9, 7e9],
        plane_k: Some(7),
        ntff: Some(fdtd::patches(sim, 2).gpu),
        fused_e,
        check_every: 250,
        ..Default::default()
    }
}

fn worst(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    let scale = a.iter().map(|v| v.abs()).fold(0.0f32, f32::max);
    assert!(scale > 0.0);
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0f32, f32::max) / scale
}

fn agree(a: &Record, b: &Record, limit: f32, what: &str) {
    assert_eq!(a.steps, b.steps, "{what}");
    for (name, x, y) in [
        ("series", &a.series, &b.series),
        ("plane", &a.plane, &b.plane),
        ("ntff", &a.ntff, &b.ntff),
    ] {
        let w = worst(x, y);
        eprintln!("{what} {name}: {w:.2e}");
        assert!(w < limit, "{what} {name}: {w}");
    }
}

#[test]
fn every_feature_steps_the_same_on_the_cpu_and_the_gpu() {
    let sim = everything();
    assert!(!sim.sheets.is_empty() && !sim.debye_edges.is_empty() && !sim.inductors.is_empty());
    let pulse = Pulse { f0: 5e9, fc: 5e9 };
    let steps = 1500;
    let cpu =
        fdtd::run_cpu(&sim, 0, &pulse, &extras(&sim, steps, false), &mut |_, _| true).unwrap();
    assert!(cpu.series.iter().all(|v| v.is_finite()));
    let _one = DEVICE.lock().unwrap_or_else(|e| e.into_inner());
    let Some(g) = gpu() else {
        return;
    };
    let apart =
        fdtd::run_blocking(&g, &sim, 0, &pulse, &extras(&sim, steps, false), &mut |_, _| true)
            .unwrap();
    agree(&cpu, &apart, 1e-5, "cpu against gpu");
    let fused =
        fdtd::run_blocking(&g, &sim, 0, &pulse, &extras(&sim, steps, true), &mut |_, _| true)
            .unwrap();
    agree(&apart, &fused, 1e-5, "separate against fused");
    let other =
        fdtd::run_blocking(&g, &sim, 1, &pulse, &extras(&sim, steps, false), &mut |_, _| true)
            .unwrap();
    let cpu_other =
        fdtd::run_cpu(&sim, 1, &pulse, &extras(&sim, steps, false), &mut |_, _| true).unwrap();
    agree(&cpu_other, &other, 1e-5, "second port");
}

#[test]
fn a_run_stops_when_progress_says_so() {
    let sim = everything();
    let pulse = Pulse { f0: 5e9, fc: 5e9 };
    let mut calls = 0;
    let out = fdtd::run_cpu(&sim, 0, &pulse, &extras(&sim, 5000, false), &mut |_, _| {
        calls += 1;
        calls < 2
    });
    assert_eq!(out.err().as_deref(), Some("cancelled"));
    assert_eq!(calls, 2);
}

#[test]
fn an_empty_cavity_rings_at_its_analytic_modes() {
    let (a, b, d) = (0.10, 0.07, 0.05);
    let mut xs: Vec<f64> = Vec::new();
    let mut t = 0.0;
    let mut s = 0.0035;
    while t < a - 1e-9 {
        xs.push(t);
        t += s;
        s = (s * 1.06).min(0.0055);
    }
    xs.push(a);
    let ys: Vec<f64> = (0..=18).map(|i| b * i as f64 / 18.0).collect();
    let zs: Vec<f64> = (0..=13).map(|i| d * i as f64 / 13.0).collect();
    let grid = Grid { x: xs, y: ys, z: zs, pml: 0 };
    let n = grid.dims();
    let mats = Materials::new(&grid);
    let port = PortDef {
        name: "probe".into(),
        columns: vec![vec![Edge { comp: 2, at: [n[0] / 3, n[1] / 3, n[2] / 3] }]],
        r: 1e7,
        reference_above: false,
    };
    let sim = Sim::new(grid, &mats, &|_, _| false, &[], vec![port], &[], Media::default());
    let extras = Extras { max_steps: 40_000, decay_db: f64::INFINITY, ..Default::default() };
    let rec =
        fdtd::run_cpu(&sim, 0, &Pulse { f0: 3.0e9, fc: 2.5e9 }, &extras, &mut |_, _| true).unwrap();
    let v: Vec<f64> = rec.series.as_chunks::<2>().0.iter().map(|p| p[0] as f64).collect();
    let spectrum = |f: f64| {
        let w = 2.0 * std::f64::consts::PI * f * sim.dt;
        let (mut re, mut im) = (0.0, 0.0);
        for (s, x) in v.iter().enumerate() {
            re += x * (w * s as f64).cos();
            im -= x * (w * s as f64).sin();
        }
        (re * re + im * im).sqrt()
    };
    for (m, nn, p) in [(1, 1, 0), (2, 1, 0), (1, 1, 1)] {
        let fa = C0 / 2.0
            * ((m as f64 / a).powi(2) + (nn as f64 / b).powi(2) + (p as f64 / d).powi(2)).sqrt();
        let peak = (0..=400)
            .map(|s| fa * (0.97 + 0.06 * s as f64 / 400.0))
            .max_by(|x, y| spectrum(*x).total_cmp(&spectrum(*y)))
            .unwrap();
        let err = (peak / fa - 1.0).abs();
        assert!(err < 0.004, "mode {m}{nn}{p}: {peak:.4e} vs {fa:.4e} ({:.3}%)", err * 100.0);
    }
}
