use sims_antenna::geometry::{Geometry, WireGeometry};
use sims_antenna::solve::{Prepared, sweep};

#[test]
fn a_sweep_reports_its_frequencies() {
    let lam = 1000.0;
    let geo = WireGeometry::new(vec![vec![[-240.0, 0.0, 0.0], [240.0, 0.0, 0.0]]], [0.0; 3]);
    let prepared = Prepared::new(&Geometry::Wire(geo), lam, 2.0, 2000);
    let seen = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let into = seen.clone();
    let listening = sims_progress::listen_here(move |r| into.lock().unwrap().push(r.stage.clone()));
    let points = sweep(&prepared, 3e8 / lam * 1e3, 0.1, 5, |_| {}, || false).expect("sweeps");
    drop(listening);
    assert_eq!(points.len(), 5);
    assert!(
        seen.lock().unwrap().iter().any(|s| s == "frequency sweep"),
        "{:?}",
        seen.lock().unwrap()
    );
}
