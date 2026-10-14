use sims_circuit::{CircuitSpec, Drive, Elem, ElemKind, Sweep, Wave, run};

fn element(name: &str, kind: ElemKind, a: Option<usize>, b: Option<usize>) -> Elem {
    Elem { name: name.into(), kind, nodes: vec![a, b] }
}

fn source(wave: Wave, dc: f64, ac: f64) -> Drive {
    Drive { net: "in".into(), p: Some(0), n: None, current: false, dc, ac, r: 0.0, wave }
}

fn rc(wave: Wave, dc: f64, ac: f64) -> CircuitSpec {
    CircuitSpec {
        schematic: "rc".into(),
        nets: vec!["in".into(), "out".into()],
        elements: vec![
            element("R1", ElemKind::Resistor { r: 1e3 }, Some(0), Some(1)),
            element("C1", ElemKind::Capacitor { c: 1e-6 }, Some(1), None),
        ],
        drives: vec![source(wave, dc, ac)],
        ac: None,
        tran: None,
        noise: None,
        checks: Vec::new(),
        record: vec![1],
        temperature: 27.0,
        netlist_hash: 0,
    }
}

#[test]
fn a_divider_halves_its_input() {
    let mut spec = rc(Wave::Dc, 5.0, 0.0);
    spec.elements[1] = element("R2", ElemKind::Resistor { r: 1e3 }, Some(1), None);
    let res = run(&spec, "divider", 0);
    assert!(res.failures.is_empty(), "{:?}", res.failures);
    assert!((res.op["out"] - 2.5).abs() < 1e-6, "{}", res.op["out"]);
    assert!((res.currents["R1"] - 2.5e-3).abs() < 1e-9, "{}", res.currents["R1"]);
}

#[test]
fn an_rc_low_pass_is_3db_down_at_its_corner() {
    let mut spec = rc(Wave::Dc, 0.0, 1.0);
    let corner = 1.0 / (std::f64::consts::TAU * 1e-3);
    spec.ac = Some(Sweep { from: corner, to: corner, per_decade: 1 });
    let res = run(&spec, "rc", 0);
    let db = res.gain_db[0].values[0];
    assert!((db + 3.0103).abs() < 1e-3, "{db} dB at {} Hz", res.freqs[0]);
    assert!((res.phase_deg[0].values[0] + 45.0).abs() < 1e-6);
}

#[test]
fn an_rc_charges_to_63_percent_in_one_time_constant() {
    let step = Wave::Pulse {
        low: 0.0,
        high: 1.0,
        delay: 0.0,
        rise: 1e-9,
        fall: 1e-9,
        width: 1.0,
        period: 2.0,
    };
    let mut spec = rc(step, 0.0, 0.0);
    spec.tran = Some((3e-3, 1e-6));
    let res = run(&spec, "rc", 0);
    assert!(res.failures.is_empty(), "{:?}", res.failures);
    let k = res.times.iter().position(|t| *t >= 1e-3).unwrap();
    let v = res.waves[0].values[k];
    assert!((v - (1.0 - (-1.0f64).exp())).abs() < 0.01, "{v} V at {} s", res.times[k]);
}

fn stages_of(run: impl FnOnce()) -> Vec<String> {
    let seen = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let into = seen.clone();
    let listening = sims_progress::listen(move |r| into.lock().unwrap().push(r.stage.clone()));
    run();
    drop(listening);
    seen.lock().unwrap().clone()
}

#[test]
fn a_run_reports_its_analyses() {
    let step = Wave::Pulse {
        low: 0.0,
        high: 1.0,
        delay: 0.0,
        rise: 1e-9,
        fall: 1e-9,
        width: 1.0,
        period: 2.0,
    };
    let mut spec = rc(step, 0.0, 1.0);
    spec.tran = Some((1e-3, 1e-5));
    spec.ac = Some(Sweep { from: 10.0, to: 1e4, per_decade: 5 });
    let stages = stages_of(|| {
        run(&spec, "rc", 0);
    });
    for want in ["finding the operating point", "AC sweep", "transient"] {
        assert!(stages.iter().any(|s| s == want), "{want} missing from {stages:?}");
    }
}

fn follower(slew: f64, wave: Wave, dc: f64) -> CircuitSpec {
    let opamp = ElemKind::Opamp {
        gain: 1e5,
        gbw: 10e6,
        rout: 1.0,
        vos: 0.0,
        en: 0.0,
        fc: 0.0,
        headroom_high: 0.0,
        headroom_low: 0.0,
        slew,
        ib: 0.0,
        psrr: 0.0,
    };
    CircuitSpec {
        schematic: "follower".into(),
        nets: vec!["in".into(), "out".into(), "vp".into()],
        elements: vec![Elem {
            name: "U1".into(),
            kind: opamp,
            nodes: vec![Some(0), Some(1), Some(1), Some(2), None],
        }],
        drives: vec![
            Drive {
                net: "in".into(),
                p: Some(0),
                n: None,
                current: false,
                dc,
                ac: 0.0,
                r: 0.0,
                wave,
            },
            Drive {
                net: "vp".into(),
                p: Some(2),
                n: None,
                current: false,
                dc: 10.0,
                ac: 0.0,
                r: 0.0,
                wave: Wave::Dc,
            },
        ],
        ac: None,
        tran: None,
        noise: None,
        checks: Vec::new(),
        record: vec![1],
        temperature: 27.0,
        netlist_hash: 0,
    }
}

fn at(res: &sims_circuit::CircuitResult, t: f64) -> f64 {
    let k = res.times.iter().position(|s| *s >= t).unwrap();
    res.waves[0].values[k]
}

#[test]
fn a_slew_limited_follower_ramps_at_its_slew_rate() {
    let step = Wave::Pulse {
        low: 2.0,
        high: 7.0,
        delay: 1e-6,
        rise: 1e-8,
        fall: 1e-8,
        width: 1.0,
        period: 2.0,
    };
    let mut spec = follower(1e6, step, 2.0);
    spec.tran = Some((10e-6, 2e-8));
    let res = run(&spec, "slew", 0);
    assert!(res.failures.is_empty(), "{:?}", res.failures);
    let rate = (at(&res, 5e-6) - at(&res, 3e-6)) / 2e-6;
    assert!((rate / 1e6 - 1.0).abs() < 0.005, "{} V/us", rate / 1e6);
}

#[test]
fn a_follower_tracks_inside_its_rails_and_clips_at_them() {
    for (input, want) in [(1.0, 1.0), (5.0, 5.0), (8.5, 8.5), (12.0, 10.0), (-2.0, 0.0)] {
        let res = run(&follower(0.0, Wave::Dc, input), "follower", 0);
        let out = res.op["out"];
        assert!((out - want).abs() < 0.02, "{input} V in, {out} V out");
    }
}

fn mains_with_third(third_rms: f64, class_b: bool) -> CircuitSpec {
    let sine = |amplitude: f64, frequency: f64| Wave::Sine {
        offset: 0.0,
        amplitude,
        frequency,
        delay: 0.0,
    };
    let drive = |net: &str, p: usize, n: Option<usize>, wave: Wave| Drive {
        net: net.into(),
        p: Some(p),
        n,
        current: false,
        dc: 0.0,
        ac: 0.0,
        r: 0.0,
        wave,
    };
    CircuitSpec {
        schematic: "mains".into(),
        nets: vec!["mid".into(), "top".into()],
        elements: vec![element("R1", ElemKind::Resistor { r: 1.0 }, Some(1), None)],
        drives: vec![
            drive("mid", 0, None, sine(325.0, 50.0)),
            drive("top", 1, Some(0), sine(third_rms * std::f64::consts::SQRT_2, 150.0)),
        ],
        ac: None,
        tran: Some((0.1, 20e-6)),
        noise: None,
        checks: vec![sims_circuit::Check {
            label: "harmonics(R1)".into(),
            measure: sims_circuit::Measure::Harmonics { element: 0, mains: 50.0, class_b },
            at: sims_circuit::At::Window(0.02, 0.1),
            min: None,
            max: Some(100.0),
        }],
        record: vec![1],
        temperature: 27.0,
        netlist_hash: 0,
    }
}

#[test]
fn harmonic_currents_read_against_the_class_a_and_b_limits() {
    let a = run(&mains_with_third(1.15, false), "a", 0);
    assert!(a.failures.is_empty(), "{:?}", a.failures);
    assert!((a.checks[0].value - 50.0).abs() < 0.5, "{:?}", a.checks[0]);
    let rms = &a.harmonics[0].values;
    assert!((rms[0] - 325.0 / std::f64::consts::SQRT_2).abs() < 0.5, "fundamental {}", rms[0]);
    assert!(rms[1].abs() < 0.01 && rms[4].abs() < 0.01, "{:?}", &rms[..6]);
    let b = run(&mains_with_third(1.15, true), "b", 0);
    assert!((b.checks[0].value - 100.0 / 3.0).abs() < 0.5, "{:?}", b.checks[0]);
    let over = run(&mains_with_third(3.0, false), "over", 0);
    assert!(!over.passed && over.checks[0].value > 100.0, "{:?}", over.checks[0]);
}

fn line(frequency: f64, dbuv: f64, average: bool) -> CircuitSpec {
    let rms = 1e-6 * 10f64.powf(dbuv / 20.0);
    let wave = Wave::Sine {
        offset: 0.0,
        amplitude: rms * std::f64::consts::SQRT_2,
        frequency,
        delay: 0.0,
    };
    CircuitSpec {
        schematic: "lisn".into(),
        nets: vec!["port".into()],
        elements: vec![element("R50", ElemKind::Resistor { r: 50.0 }, Some(0), None)],
        drives: vec![Drive {
            net: "port".into(),
            p: Some(0),
            n: None,
            current: false,
            dc: 0.0,
            ac: 0.0,
            r: 0.0,
            wave,
        }],
        ac: None,
        tran: Some((500e-6, 10e-9)),
        noise: None,
        checks: vec![sims_circuit::Check {
            label: "emission(port)".into(),
            measure: sims_circuit::Measure::Emission { node: Some(0), class_a: false, average },
            at: sims_circuit::At::Window(0.0, 500e-6),
            min: None,
            max: Some(0.0),
        }],
        record: vec![0],
        temperature: 27.0,
        netlist_hash: 0,
    }
}

#[test]
fn a_line_reads_its_level_against_the_cispr_11_class_b_limits() {
    let one = run(&line(1e6, 60.0, false), "1MHz", 0);
    assert!(!one.passed && one.emission_freqs.len() > 6000, "{}", one.emission_freqs.len());
    assert!((one.checks[0].value - 4.0).abs() < 0.3, "{:?}", one.checks[0]);
    let average = run(&line(1e6, 60.0, true), "1MHz average", 0);
    assert!((average.checks[0].value - 14.0).abs() < 0.3, "{:?}", average.checks[0]);
    let low = run(&line(200e3, 70.0, false), "200kHz", 0);
    let limit = 66.0 - 10.0 * (200f64 / 150.0).log10() / (500f64 / 150.0).log10();
    assert!(
        (low.checks[0].value - (70.0 - limit)).abs() < 0.3,
        "{:?} against {limit}",
        low.checks[0]
    );
    let quiet = run(&line(1e6, 40.0, false), "quiet", 0);
    assert!(quiet.passed, "{:?}", quiet.checks[0]);
}
