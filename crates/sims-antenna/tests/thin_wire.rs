use sims_antenna::fdtd::{Model, Options, Port, Wire, run};
use sims_antenna::geometry::WireGeometry;
use sims_antenna::solve::{model_for, solve_at};

const LEN: f64 = 0.141;
const GAP: f64 = 0.003;

fn dipole(radius: f64) -> Model {
    let half = |s: f64| Wire { a: [0.0, 0.0, s * GAP / 2.0], b: [0.0, 0.0, s * LEN / 2.0], radius };
    Model {
        wires: vec![half(-1.0), half(1.0)],
        port: Some(Port { a: [0.0, 0.0, -GAP / 2.0], b: [0.0, 0.0, GAP / 2.0], r: 50.0 }),
        f0: 1.0e9,
        span: 0.2,
        fine: GAP,
        ..Default::default()
    }
}

fn fdtd_resonance(radius: f64) -> f64 {
    let out = run(&dipole(radius), &Options::default(), |_, _| true).unwrap();
    let pts: Vec<_> = out.sweep.freqs.iter().zip(&out.sweep.z).collect();
    let w = pts.windows(2).find(|w| w[0].1.im < 0.0 && w[1].1.im >= 0.0).expect("a resonance");
    let (f0, f1, x0, x1) = (*w[0].0, *w[1].0, w[0].1.im, w[1].1.im);
    f0 - x0 * (f1 - f0) / (x1 - x0)
}

fn mom_resonance(radius: f64) -> f64 {
    let x = |f: f64| {
        let lam = 299_792_458.0 / f * 1e3;
        let ends = vec![[0.0, 0.0, -LEN * 500.0], [0.0, 0.0, LEN * 500.0]];
        let g = WireGeometry::new(vec![ends], [0.0; 3]);
        solve_at(&model_for(&g, lam, 2.0 * radius * 1e3, 400), lam, false).z.im
    };
    let (mut lo, mut hi) = (0.8e9, 1.2e9);
    for _ in 0..30 {
        let m = (lo + hi) / 2.0;
        if x(m) < 0.0 { lo = m } else { hi = m }
    }
    lo
}

#[test]
fn a_thicker_wire_resonates_lower_by_what_the_moment_method_says() {
    let radii = [0.1e-3, 1.2e-3];
    let fdtd = radii.map(fdtd_resonance);
    let mom = radii.map(mom_resonance);
    for k in 0..2 {
        assert!((fdtd[k] / mom[k] - 1.0).abs() < 0.025, "{} against {}", fdtd[k], mom[k]);
    }
    let ratio = (fdtd[0] - fdtd[1]) / (mom[0] - mom[1]);
    assert!((0.6..1.2).contains(&ratio), "{fdtd:?} against {mom:?}");
}

#[test]
fn a_wire_must_run_along_an_axis_and_fit_its_cell() {
    let mut slanted = dipole(0.1e-3);
    slanted.wires[1].b[0] = 0.01;
    let err = |m: &Model| run(m, &Options::default(), |_, _| true).err().unwrap_or_default();
    assert!(err(&slanted).contains("along one axis"), "{}", err(&slanted));
    for radius in [2.0e-3, 0.01e-3] {
        let e = err(&dipole(radius));
        assert!(e.contains("from a 55th to half of its cell"), "{e}");
    }
}
