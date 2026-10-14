use sims_logic::{
    Cell, Circuit, GateOp, Input, Level, OnViolation, Output, Prim, Stimulus, Wave, simulate,
};

#[test]
fn a_run_reports_simulated_time() {
    let inverter = Cell {
        part: "U1".into(),
        prim: Prim::Gate { op: GateOp::And, invert: true },
        input_names: vec!["A".into()],
        inputs: vec![Input::Net { net: 0, invert: false }],
        outputs: vec![Output { net: Some(1), invert: false, delay: 1000 }],
        weak: false,
        open_drain: false,
        setup: 0,
        hold: 0,
        recovery: 0,
        removal: 0,
    };
    let circuit = Circuit { nets: vec!["A".into(), "Y".into()], cells: vec![inverter] };
    let clock = Stimulus {
        net: 0,
        name: "A".into(),
        wave: Wave::Clock { period: 20_000, high: 10_000, phase: 0 },
    };
    let seen = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let into = seen.clone();
    let listening = sims_progress::listen_here(move |r| into.lock().unwrap().push(r.stage.clone()));
    let run = simulate(&circuit, &[clock], 1_000_000, OnViolation::X);
    drop(listening);
    assert_eq!(run.changes[1].last().map(|c| c.1), Some(Level::H));
    assert!(seen.lock().unwrap().iter().any(|s| s == "logic simulation"));
}
