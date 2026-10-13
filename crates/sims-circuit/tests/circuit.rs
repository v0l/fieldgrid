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
