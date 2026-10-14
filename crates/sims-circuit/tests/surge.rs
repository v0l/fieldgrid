use sims_circuit::{CircuitSpec, Drive, Elem, ElemKind, Surge, Wave, run};

fn crossing(wave: &Wave, level: f64, from: f64, to: f64, rising: bool) -> f64 {
    let n = 200_000;
    let at = |k: usize| from + (to - from) * k as f64 / n as f64;
    (0..n)
        .map(at)
        .find(|&t| {
            let v = wave.at(0.0, t);
            if rising { v >= level } else { v <= level }
        })
        .unwrap_or(f64::NAN)
}

fn peak_time(wave: &Wave, to: f64) -> (f64, f64) {
    let n = 200_000;
    (0..n)
        .map(|k| to * k as f64 / n as f64)
        .map(|t| (t, wave.at(0.0, t)))
        .fold((0.0, f64::MIN), |a, b| if b.1 > a.1 { b } else { a })
}

#[test]
fn the_open_circuit_surge_is_1_2_over_50() {
    let wave = Wave::Surge { offset: 0.0, peak: 1.0, delay: 0.0, shape: Surge::Open };
    let (t, peak) = peak_time(&wave, 10e-6);
    assert!((peak - 1.0).abs() < 1e-6, "{peak}");
    let front = 1.67 * (crossing(&wave, 0.9, 0.0, t, true) - crossing(&wave, 0.3, 0.0, t, true));
    let width = crossing(&wave, 0.5, t, 200e-6, false) - crossing(&wave, 0.5, 0.0, t, true);
    assert!((front / 1.2e-6 - 1.0).abs() < 0.05, "front {front}");
    assert!((width / 50e-6 - 1.0).abs() < 0.05, "width {width}");
}

#[test]
fn the_short_circuit_surge_is_8_over_20() {
    let wave = Wave::Surge { offset: 0.0, peak: 1.0, delay: 0.0, shape: Surge::Short };
    let (t, peak) = peak_time(&wave, 40e-6);
    assert!((peak - 1.0).abs() < 1e-6, "{peak}");
    let front = 1.25 * (crossing(&wave, 0.9, 0.0, t, true) - crossing(&wave, 0.1, 0.0, t, true));
    let duration =
        1.18 * (crossing(&wave, 0.5, t, 200e-6, false) - crossing(&wave, 0.5, 0.0, t, true));
    assert!((front / 8e-6 - 1.0).abs() < 0.05, "front {front}");
    assert!((duration / 20e-6 - 1.0).abs() < 0.05, "duration {duration}");
}

#[test]
fn a_burst_is_a_train_of_5_over_50_ns_pulses() {
    let wave = Wave::Burst {
        offset: 0.0,
        peak: 1.0,
        delay: 0.0,
        rate: 5e3,
        length: 15e-3,
        period: 300e-3,
    };
    let pulse =
        Wave::Burst { offset: 0.0, peak: 1.0, delay: 0.0, rate: 5e3, length: 1e-4, period: 0.0 };
    let (t, peak) = peak_time(&pulse, 100e-9);
    assert!((peak - 1.0).abs() < 1e-3, "{peak}");
    let rise = crossing(&pulse, 0.9, 0.0, t, true) - crossing(&pulse, 0.1, 0.0, t, true);
    let width = crossing(&pulse, 0.5, t, 400e-9, false) - crossing(&pulse, 0.5, 0.0, t, true);
    assert!((rise - 5e-9).abs() < 1.5e-9, "rise {rise}");
    assert!((width - 50e-9).abs() < 15e-9, "width {width}");
    let starts = wave
        .breaks(0.6)
        .iter()
        .filter(|b| {
            let since = (*b % 200e-6 + 200e-6) % 200e-6;
            since < 1e-12 || 200e-6 - since < 1e-12
        })
        .count();
    assert_eq!(starts, 2 * 75 + 1, "pulse starts in two bursts");
    assert!(wave.at(0.0, 20e-3 + t).abs() < 1e-12, "quiet between bursts");
    assert!((wave.at(0.0, 300e-3 + 10e-3 + t) - 1.0).abs() < 1e-3, "second burst");
}

fn mov() -> ElemKind {
    ElemKind::Varistor { v1ma: 430.0, alpha: 30.0, c: 0.0 }
}

fn spec(drive: Drive) -> CircuitSpec {
    CircuitSpec {
        schematic: "mov".into(),
        nets: vec!["a".into()],
        elements: vec![Elem { name: "RV1".into(), kind: mov(), nodes: vec![Some(0), None] }],
        drives: vec![drive],
        ac: None,
        tran: None,
        noise: None,
        checks: Vec::new(),
        record: vec![0],
        temperature: 27.0,
        netlist_hash: 0,
    }
}

#[test]
fn a_varistor_follows_its_power_law() {
    for amps in [1e-3, 1.0, 100.0, -10.0] {
        let drive = Drive {
            net: "a".into(),
            p: Some(0),
            n: None,
            current: true,
            dc: amps,
            ac: 0.0,
            r: 0.0,
            wave: Wave::Dc,
        };
        let res = run(&spec(drive), "mov", 0);
        assert!(res.failures.is_empty(), "{:?}", res.failures);
        let want = 430.0 * (amps.abs() / 1e-3).powf(1.0 / 30.0) * amps.signum();
        let got = res.op["a"];
        assert!((got / want - 1.0).abs() < 1e-6, "{amps} A: {got} V against {want}");
    }
}

#[test]
fn a_surge_through_two_ohms_is_clamped_by_a_varistor() {
    let peak = 2000.0;
    let drive = Drive {
        net: "a".into(),
        p: Some(0),
        n: None,
        current: false,
        dc: 0.0,
        ac: 0.0,
        r: 2.0,
        wave: Wave::Surge { offset: 0.0, peak, delay: 1e-6, shape: Surge::Open },
    };
    let mut circuit = spec(drive);
    circuit.tran = Some((100e-6, 1e-6));
    let res = run(&circuit, "surge", 0);
    assert!(res.failures.is_empty(), "{:?}", res.failures);
    let clamp = res.waves[0].values.iter().copied().fold(f64::MIN, f64::max);
    let current = |v: f64| 1e-3 * (v / 430.0).powf(30.0);
    let (mut lo, mut hi) = (430.0, peak);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if current(mid) * 2.0 + mid > peak { hi = mid } else { lo = mid }
    }
    assert!((clamp / lo - 1.0).abs() < 0.01, "clamped at {clamp} V against {lo} V");
}
