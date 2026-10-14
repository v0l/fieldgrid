use sims_acoustics::{Air, BANDS, Rays, Response, Room, shoebox, trace};

fn close(got: f64, want: f64, tol: f64) {
    assert!((got / want - 1.0).abs() < tol, "{got} against {want}");
}

#[test]
fn air_absorption_matches_iso_9613_1_table_1() {
    let air = |humidity| Air { celsius: -20.0, humidity, kilopascal: 101.325 };
    let exact = |k: i32| 1000.0 * 10f64.powf(k as f64 / 10.0);
    for (humidity, k, per_km) in [
        (50.0, 0, 9.14),
        (50.0, 6, 13.2),
        (50.0, -6, 2.46),
        (10.0, 0, 1.65),
        (10.0, 9, 10.9),
        (100.0, 0, 11.1),
        (30.0, -3, 4.38),
    ] {
        close(air(humidity).decibels_per_metre(exact(k)) * 1000.0, per_km, 0.006);
    }
}

fn exponential(t60: f64, bin: f64) -> Response {
    let decay = 6.0 * std::f64::consts::LN_10 / t60;
    let curve: Vec<f64> = (0..(4.0 * t60 / bin) as usize)
        .map(|i| (-decay * i as f64 * bin).exp() - (-decay * (i + 1) as f64 * bin).exp())
        .collect();
    Response { bin, direct: 0.0, energy: std::array::from_fn(|_| curve.clone()) }
}

#[test]
fn an_exponential_decay_reads_its_own_figures() {
    let t = 1.2;
    let r = exponential(t, 1e-4);
    close(r.t30(3).unwrap(), t, 1e-3);
    close(r.edt(3).unwrap(), t, 1e-3);
    let k = 6.0 * std::f64::consts::LN_10 / t;
    close(r.definition(3), 1.0 - (-k * 0.05).exp(), 1e-3);
    close(r.clarity(3, 80.0), 10.0 * ((k * 0.08).exp() - 1.0).log10(), 1e-3);
    for f in [0.63, 2.0, 12.5] {
        let m = 1.0 / (1.0 + (2.0 * std::f64::consts::PI * f * t / 13.8).powi(2)).sqrt();
        close(r.modulation(3, f), m, 2e-3);
    }
    assert!(exponential(0.3, 1e-4).sti() > exponential(1.2, 1e-4).sti());
    let mut dry = exponential(1.0, 1e-3);
    dry.energy.iter_mut().for_each(|e| *e = vec![1.0]);
    close(dry.sti(), 1.0, 1e-9);
}

fn room(absorption: f64) -> Room {
    Room {
        surfaces: shoebox([10.0, 7.0, 4.0], [absorption; BANDS.len()], 1.0),
        air: Air::default(),
        air_absorption: false,
    }
}

#[test]
fn a_diffuse_shoebox_decays_between_eyring_and_sabine() {
    for (absorption, near_eyring) in [(0.05, 0.015), (0.1, 0.025), (0.3, 0.12)] {
        let r = room(absorption);
        close(r.volume(), 280.0, 1e-12);
        close(r.area(), 276.0, 1e-12);
        let rays = Rays { count: 4000, seconds: 1.5 * r.eyring(3), ..Default::default() };
        let got = trace(&r, [2.1, 3.3, 1.7], &[[7.4, 2.2, 1.2]], &rays).unwrap();
        let t = got[0].t30(3).unwrap();
        assert!(
            r.eyring(3) < t && t < r.sabine(3),
            "{t} outside {} to {}",
            r.eyring(3),
            r.sabine(3)
        );
        close(t, r.eyring(3), near_eyring);
    }
}

#[test]
fn air_absorption_shortens_the_treble() {
    let mut r = room(0.1);
    r.air_absorption = true;
    let rays = Rays { count: 3000, seconds: 2.5, ..Default::default() };
    let got = &trace(&r, [2.1, 3.3, 1.7], &[[7.4, 2.2, 1.2]], &rays).unwrap()[0];
    let (low, high) = (got.t30(0).unwrap(), got.t30(6).unwrap());
    close(high, r.eyring(6), 0.05);
    assert!(high < 0.8 * low, "{high} against {low}");
}

#[test]
fn tracing_reports_its_rays() {
    let seen = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let into = seen.clone();
    let listening = sims_progress::listen_here(move |r| into.lock().unwrap().push(r.stage.clone()));
    let rays = Rays { count: 50, seconds: 0.2, ..Default::default() };
    trace(&room(0.3), [2.0, 2.0, 2.0], &[[5.0, 3.0, 2.0]], &rays).unwrap();
    drop(listening);
    assert!(seen.lock().unwrap().iter().any(|s| s == "tracing rays"));
}
