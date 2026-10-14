use sims_thermal::{Block, Boundary, Grid, Model, lines};

fn plane(axis: usize, at: f64, environment: usize, resistance: f64) -> Boundary {
    let (mut low, mut high) = ([-1e9; 3], [1e9; 3]);
    low[axis] = at;
    high[axis] = at;
    Boundary { low, high, environment, resistance }
}

#[test]
fn case_1_half_square_column() {
    let x = lines(&[0.0, 0.25, 0.5, 0.75, 1.0], 0.0125, 1.0);
    let y = lines(&[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0], 0.0125, 1.0);
    let model = Model {
        grid: Grid::plane(x, y),
        blocks: vec![Block { low: [0.0; 3], high: [1.0, 2.0, 1.0], conductivity: 1.0 }],
        boundaries: vec![plane(1, 2.0, 0, 0.0), plane(0, 0.0, 1, 0.0), plane(1, 0.0, 1, 0.0)],
    };
    let prepared = model.prepare().unwrap();
    let s = prepared.solve(&[20.0, 0.0]).unwrap();
    let iso = [
        [9.7, 13.4, 14.7, 15.1],
        [5.3, 8.6, 10.3, 10.8],
        [3.2, 5.6, 7.0, 7.5],
        [2.0, 3.6, 4.7, 5.0],
        [1.3, 2.3, 3.0, 3.2],
        [0.7, 1.4, 1.8, 1.9],
        [0.3, 0.6, 0.8, 0.9],
    ];
    let analytic = |x: f64, y: f64| {
        (1..400)
            .step_by(2)
            .map(|n| {
                let k = n as f64 * std::f64::consts::PI / 2.0;
                80.0 / (n as f64 * std::f64::consts::PI)
                    * (k * x).sin()
                    * (k * (y - 2.0)).exp()
                    * (1.0 - (-2.0 * k * y).exp())
                    / (1.0 - (-4.0 * k).exp())
            })
            .sum::<f64>()
    };
    for (r, row) in iso.iter().enumerate() {
        for (c, want) in row.iter().enumerate() {
            let (x, y) = (0.25 * (c + 1) as f64, 1.75 - 0.25 * r as f64);
            let got = prepared.temperature(&s, [x, y, 0.5]).unwrap();
            assert!((got - want).abs() <= 0.1, "row {r} column {c}: {got} against {want}");
            assert!(
                (got - analytic(x, y)).abs() < 0.02,
                "row {r} column {c}: {got} against {}",
                analytic(x, y)
            );
        }
    }
}

#[test]
fn case_2_composite_wall() {
    let mm = 1e-3;
    let x = lines(&[0.0, 1.5 * mm, 15.0 * mm, 500.0 * mm], 0.25 * mm, 1.15);
    let y = lines(&[0.0, 1.5 * mm, 35.0 * mm, 36.5 * mm, 41.5 * mm, 47.5 * mm], 0.25 * mm, 1.15);
    let block = |x0: f64, y0: f64, x1: f64, y1: f64, conductivity| Block {
        low: [x0 * mm, y0 * mm, 0.0],
        high: [x1 * mm, y1 * mm, 1.0],
        conductivity,
    };
    let model = Model {
        grid: Grid::plane(x, y),
        blocks: vec![
            block(0.0, 0.0, 500.0, 47.5, 0.029),
            block(0.0, 41.5, 500.0, 47.5, 1.15),
            block(0.0, 36.5, 15.0, 41.5, 0.12),
            block(0.0, 0.0, 500.0, 1.5, 230.0),
            block(0.0, 0.0, 1.5, 36.5, 230.0),
            block(0.0, 35.0, 15.0, 36.5, 230.0),
        ],
        boundaries: vec![plane(1, 47.5 * mm, 0, 0.06), plane(1, 0.0, 1, 0.11)],
    };
    let prepared = model.prepare().unwrap();
    let s = prepared.solve(&[0.0, 20.0]).unwrap();
    assert!((s.flows[1] - 9.5).abs() <= 0.1, "{} W/m", s.flows[1]);
    let iterated = prepared.solve_iterative(&[0.0, 20.0], None, 1e-10).unwrap();
    assert!((iterated.flows[1] / s.flows[1] - 1.0).abs() < 1e-6, "{} W/m", iterated.flows[1]);
    if let Ok(gpu) = sims_field::Gpu::open_blocking(&Default::default()) {
        let g = prepared.solve_iterative(&[0.0, 20.0], Some(&gpu), 1e-6).unwrap();
        assert!((g.flows[1] / s.flows[1] - 1.0).abs() < 2e-3, "{} W/m on {}", g.flows[1], gpu.name);
    }
    assert!((s.flows[0] + s.flows[1]).abs() < 1e-6 * s.flows[1]);
    for (name, x, y, want) in [
        ("A", 0.0, 47.5, 7.1),
        ("B", 500.0, 47.5, 0.8),
        ("C", 0.0, 41.5, 7.9),
        ("D", 15.0, 41.5, 6.3),
        ("E", 500.0, 41.5, 0.8),
        ("F", 0.0, 36.5, 16.4),
        ("G", 15.0, 36.5, 16.3),
        ("H", 0.0, 0.0, 16.8),
        ("I", 500.0, 0.0, 18.3),
    ] {
        let got = prepared.temperature(&s, [x * mm, y * mm, 0.5]).unwrap();
        assert!((got - want).abs() <= 0.1, "{name}: {got} against {want}");
    }
}

#[test]
fn grid_lines_keep_their_stops_and_grow_evenly() {
    for (stops, fine, grow) in
        [(vec![0.0, 0.25, 1.0], 0.0125, 1.0), (vec![0.0, 0.0015, 0.5], 0.00025, 1.15)]
    {
        let l = lines(&stops, fine, grow);
        assert!(l.windows(2).all(|w| w[1] - w[0] > 0.4 * fine), "{l:?}");
        assert!(stops.iter().all(|s| l.iter().any(|v| (v - s).abs() < 1e-15)));
        let widest = l.windows(2).map(|w| w[1] - w[0]).fold(0.0, f64::max);
        assert!(widest < (stops[2] - stops[1]) / 2.0, "{widest}");
    }
}

#[test]
fn coupling_coefficients_of_a_slab_add_its_resistances() {
    let model = Model {
        grid: Grid::plane(lines(&[0.0, 1.0], 0.1, 1.0), lines(&[0.0, 0.2], 0.01, 1.0)),
        blocks: vec![Block { low: [0.0; 3], high: [1.0, 0.2, 1.0], conductivity: 0.04 }],
        boundaries: vec![plane(1, 0.0, 0, 0.13), plane(1, 0.2, 1, 0.04)],
    };
    let l = model.prepare().unwrap().coupling().unwrap();
    let u = 1.0 / (0.13 + 0.2 / 0.04 + 0.04);
    for (i, j, want) in [(0, 0, u), (1, 1, u), (0, 1, u), (1, 0, u)] {
        assert!((l[i][j] / want - 1.0).abs() < 1e-6, "{l:?}");
    }
}

#[test]
fn solving_reports_its_stages() {
    let seen = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let into = seen.clone();
    let listening = sims_progress::listen_here(move |r| into.lock().unwrap().push(r.stage.clone()));
    let model = Model {
        grid: Grid::plane(vec![0.0, 1.0], vec![0.0, 1.0]),
        blocks: vec![Block { low: [0.0; 3], high: [1.0; 3], conductivity: 1.0 }],
        boundaries: vec![plane(0, 0.0, 0, 0.1), plane(0, 1.0, 1, 0.1)],
    };
    model.prepare().unwrap().solve(&[0.0, 1.0]).unwrap();
    drop(listening);
    let seen = seen.lock().unwrap();
    assert!(seen.iter().any(|s| s == "solving the conduction grid"), "{seen:?}");
}

fn case_4(fine: f64) -> (f64, f64) {
    let x = lines(&[0.0, 0.45, 0.55, 1.0], fine, 1.2);
    let y = lines(&[0.0, 0.2, 0.6], fine, 1.2);
    let z = lines(&[0.0, 0.475, 0.525, 1.0], fine, 1.2);
    let model = Model {
        grid: Grid { x: x.clone(), y, z: z.clone() },
        blocks: vec![
            Block { low: [0.0; 3], high: [1.0, 0.2, 1.0], conductivity: 0.1 },
            Block { low: [0.45, 0.0, 0.475], high: [0.55, 0.6, 0.525], conductivity: 50.0 },
        ],
        boundaries: vec![
            plane(1, 0.0, 0, 0.1),
            plane(1, 0.2, 1, 0.1),
            Boundary {
                low: [0.45, 0.2, 0.475],
                high: [0.55, 0.6, 0.525],
                environment: 1,
                resistance: 0.1,
            },
        ],
    };
    let prepared = model.prepare().unwrap();
    let s = prepared.solve(&[0.0, 1.0]).unwrap();
    let hottest = x
        .iter()
        .flat_map(|a| z.iter().map(move |c| (*a, *c)))
        .filter_map(|(a, c)| prepared.temperature(&s, [a, 0.0, c]))
        .fold(f64::MIN, f64::max);
    (s.flows[1], hottest)
}

#[test]
fn case_4_iron_bar_through_insulation() {
    let (flow, hottest) = case_4(0.0025);
    assert!((flow / 0.540 - 1.0).abs() <= 0.01, "{flow} W");
    assert!((hottest - 0.805).abs() <= 0.005, "{hottest}");
}

fn case_3(fine: f64) -> (Vec<f64>, f64, f64, Vec<Vec<f64>>) {
    let x = lines(&[-0.1, 0.0, 0.05, 0.2, 1.2], fine, 1.25);
    let y = lines(&[-0.7, -0.1, 0.0, 0.05, 0.2, 1.2], fine, 1.25);
    let z = lines(&[0.0, 1.0, 1.15, 1.2, 2.15], fine, 1.25);
    let block = |low: [f64; 3], high: [f64; 3], conductivity| Block { low, high, conductivity };
    let wall = |low: [f64; 2], high: [f64; 2], conductivity| {
        [
            block([low[0], low[1], 0.0], [high[0], high[1], 2.15], conductivity),
            block([low[1], low[0], 0.0], [high[1], high[0], 2.15], conductivity),
        ]
    };
    let mut blocks = Vec::new();
    blocks.extend(wall([-0.1, -0.1], [1.2, 0.0], 1.0));
    blocks.extend(wall([0.0, 0.0], [1.2, 0.05], 0.04));
    blocks.extend(wall([0.05, 0.05], [1.2, 0.2], 0.7));
    blocks.push(block([0.05, -0.7, 1.0], [1.2, 1.2, 1.15], 2.5));
    blocks.push(block([0.2, 0.2, 1.15], [1.2, 1.2, 1.2], 1.0));
    let region = |low: [f64; 3], high: [f64; 3], environment, resistance| Boundary {
        low,
        high,
        environment,
        resistance,
    };
    let adiabatic = f64::INFINITY;
    let boundaries = vec![
        plane(0, -0.1, 0, 0.05),
        plane(1, -0.1, 0, 0.05),
        region([0.05, -0.7, 1.0], [1.2, -0.1, 1.15], 0, 0.05),
        region([0.2, 0.2, 0.0], [1.2, 1.2, 1.0], 1, 0.2),
        region([0.2, 0.2, 1.15], [1.2, 1.2, 2.15], 2, 0.2),
        plane(0, 1.2, 0, adiabatic),
        plane(1, 1.2, 0, adiabatic),
        plane(2, 0.0, 0, adiabatic),
        plane(2, 2.15, 0, adiabatic),
    ];
    let model = Model { grid: Grid { x, y, z }, blocks, boundaries };
    let prepared = model.prepare().unwrap();
    let s = prepared.solve(&[0.0, 20.0, 15.0]).unwrap();
    let alpha = prepared.temperature(&s, [0.2, 0.2, 1.0]).unwrap();
    let beta = prepared.temperature(&s, [0.2, 0.2, 1.2]).unwrap();
    (s.flows, alpha, beta, prepared.coupling().unwrap())
}

#[test]
fn case_3_corner_with_a_floor_slab() {
    let (flows, alpha, beta, l) = case_3(0.01);
    for (got, want) in [(-flows[0], 59.98), (flows[1], 46.09), (flows[2], 13.89)] {
        assert!((got / want - 1.0).abs() <= 0.01, "{got} W against {want}");
    }
    assert!((alpha - 11.32).abs() <= 0.1, "{alpha}");
    assert!((beta - 11.11).abs() <= 0.1, "{beta}");
    for (i, j, want) in [(1, 0, 1.781), (2, 0, 1.624), (1, 2, 2.094)] {
        assert!((l[i][j] / want - 1.0).abs() <= 0.01, "L{i}{j} {} against {want}", l[i][j]);
    }
}
