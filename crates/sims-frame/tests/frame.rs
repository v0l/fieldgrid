use sims_frame::{Case, FIXED, Material, Member, Model, PINNED, Section, buckling, modes, solve};
use std::f64::consts::PI;

const STEEL: Material = Material { e: 200e9, g: 200e9 / 2.6, density: 7850.0 };

fn member(a: usize, b: usize, section: Section) -> Member {
    Member { a, b, section, material: STEEL, up: [0.0, 0.0, 1.0], pinned: [false; 2] }
}

fn line(from: [f64; 3], to: [f64; 3], count: usize, section: Section) -> Model {
    let nodes = (0..=count)
        .map(|i| {
            let t = i as f64 / count as f64;
            std::array::from_fn(|d| from[d] + (to[d] - from[d]) * t)
        })
        .collect();
    let members = (0..count).map(|i| member(i, i + 1, section)).collect();
    Model { nodes, members, ..Default::default() }
}

fn close(got: f64, want: f64, tol: f64) {
    assert!((got / want - 1.0).abs() < tol, "{got} against {want}");
}

#[test]
fn a_cantilever_bends_and_twists_as_the_formulas_say() {
    let s = Section::rectangle(0.1, 0.2);
    let length = 3.0;
    let mut beam = line([0.0; 3], [length, 0.0, 0.0], 6, s);
    beam.supports = vec![(0, FIXED)];
    let tip = 6;
    let (p, t) = (1000.0, 200.0);
    let case = Case { nodal: vec![(tip, [0.0, -p, -p, t, 0.0, 0.0])], ..Default::default() };
    let r = solve(&beam, &case).unwrap();
    let d = r.displacement[tip];
    close(-d[2], p * length.powi(3) / (3.0 * STEEL.e * s.iz), 1e-9);
    close(-d[1], p * length.powi(3) / (3.0 * STEEL.e * s.iy), 1e-9);
    close(d[3], t * length / (STEEL.g * s.j), 1e-9);
    close(r.reaction[0][2], p, 1e-9);
    close(r.reaction[0][4], -p * length, 1e-9);
}

#[test]
fn a_simply_supported_beam_under_a_uniform_load() {
    let s = Section::rectangle(0.1, 0.3);
    let (length, w) = (6.0, 5000.0);
    let mut beam = line([0.0; 3], [length, 0.0, 0.0], 8, s);
    beam.supports = vec![
        (0, [true, true, true, true, false, false]),
        (8, [false, true, true, false, false, false]),
    ];
    let case =
        Case { uniform: (0..8).map(|i| (i, [0.0, 0.0, -w])).collect(), ..Default::default() };
    let r = solve(&beam, &case).unwrap();
    close(-r.displacement[4][2], 5.0 * w * length.powi(4) / (384.0 * STEEL.e * s.iz), 1e-9);
    close(r.reaction[0][2], w * length / 2.0, 1e-9);
    close(r.ends[3][11].abs(), w * length * length / 8.0, 1e-9);
}

#[test]
fn pinned_ends_turn_a_fixed_beam_into_a_simple_one() {
    let s = Section::rectangle(0.1, 0.3);
    let (length, w) = (6.0, 5000.0);
    let mut beam = line([0.0; 3], [length, 0.0, 0.0], 2, s);
    beam.members[0].pinned = [true, false];
    beam.members[1].pinned = [false, true];
    beam.supports = vec![(0, FIXED), (2, FIXED)];
    let case =
        Case { uniform: vec![(0, [0.0, 0.0, -w]), (1, [0.0, 0.0, -w])], ..Default::default() };
    let r = solve(&beam, &case).unwrap();
    close(-r.displacement[1][2], 5.0 * w * length.powi(4) / (384.0 * STEEL.e * s.iz), 1e-9);
    assert!(r.reaction[0][4].abs() < 1e-6 * w * length * length, "{:?}", r.reaction[0]);
}

#[test]
fn a_fixed_portal_sways_as_slope_deflection_says() {
    let (h, span, push) = (3.5, 6.0, 10_000.0);
    let column = Section { area: 1.0, ..Section::rectangle(0.3, 0.3) };
    let girder = Section { area: 1.0, ..Section::rectangle(0.3, 0.5) };
    let nodes = vec![[0.0, 0.0, 0.0], [0.0, 0.0, h], [span, 0.0, h], [span, 0.0, 0.0]];
    let side = [1.0, 0.0, 0.0];
    let members = vec![
        Member { up: side, ..member(0, 1, column) },
        member(1, 2, girder),
        Member { up: side, ..member(3, 2, column) },
    ];
    let mut portal =
        Model { nodes, members, supports: vec![(0, FIXED), (3, FIXED)], ..Default::default() };
    let case = Case { nodal: vec![(1, [push, 0.0, 0.0, 0.0, 0.0, 0.0])], ..Default::default() };
    let sway = solve(&portal, &case).unwrap().displacement[1][0];
    let (ic, ib) = (column.iz, girder.iz);
    let k = (ib / span) / (ic / h);
    let want = push * h.powi(3) / (12.0 * STEEL.e * ic * (2.0 - 6.0 / (4.0 + 6.0 * k)));
    close(sway, want, 2e-3);
    portal.supports = vec![(0, PINNED), (3, PINNED)];
    assert!(solve(&portal, &case).unwrap().displacement[1][0] > 2.0 * sway);
}

#[test]
fn columns_buckle_at_eulers_load() {
    let s = Section::tube(0.1, 0.005);
    let length = 4.0;
    let mut column = line([0.0; 3], [0.0, 0.0, length], 10, s);
    column.members.iter_mut().for_each(|m| m.up = [1.0, 0.0, 0.0]);
    let down = Case { nodal: vec![(10, [0.0, 0.0, -1.0, 0.0, 0.0, 0.0])], ..Default::default() };
    let euler = PI * PI * STEEL.e * s.iy / (length * length);
    column.supports = vec![
        (0, [true, true, true, false, false, true]),
        (10, [true, true, false, false, false, false]),
    ];
    close(buckling(&column, &down, 1).unwrap()[0].value, euler, 1e-4);
    column.supports = vec![(0, FIXED)];
    close(buckling(&column, &down, 1).unwrap()[0].value, euler / 4.0, 1e-4);
}

#[test]
fn a_cantilever_rings_at_its_beam_frequencies() {
    let s = Section::rectangle(0.05, 0.1);
    let length = 2.0;
    let mut beam = line([0.0; 3], [length, 0.0, 0.0], 20, s);
    beam.supports = vec![(0, FIXED)];
    let found = modes(&beam, 8).unwrap();
    let omega = |beta: f64, i: f64| {
        beta * beta * (STEEL.e * i / (STEEL.density * s.area * length.powi(4))).sqrt()
    };
    let weak = found[0].frequency() * 2.0 * PI;
    let strong: Vec<f64> = found.iter().map(|m| m.frequency() * 2.0 * PI).collect();
    close(weak, omega(1.87510, s.iy), 1e-4);
    assert!(strong.iter().any(|w| (w / omega(1.87510, s.iz) - 1.0).abs() < 1e-4), "{strong:?}");
    assert!(strong.iter().any(|w| (w / omega(4.69409, s.iy) - 1.0).abs() < 1e-4), "{strong:?}");
}

#[test]
fn combinations_add_like_the_loads_do() {
    let mut beam = line([0.0; 3], [3.0, 0.0, 0.0], 3, Section::rectangle(0.1, 0.2));
    beam.supports = vec![(0, FIXED)];
    let dead = Case { gravity: Some([0.0, 0.0, -9.81]), ..Default::default() };
    let live =
        Case { nodal: vec![(3, [0.0, 500.0, -2000.0, 0.0, 0.0, 0.0])], ..Default::default() };
    let both = Case { gravity: dead.gravity, nodal: live.nodal.clone(), ..Default::default() };
    let (d, l, b) =
        (solve(&beam, &dead).unwrap(), solve(&beam, &live).unwrap(), solve(&beam, &both).unwrap());
    let c = sims_frame::Static::combine(&[(&d, 1.0), (&l, 1.0)]).unwrap();
    for (x, y) in c.displacement.iter().flatten().zip(b.displacement.iter().flatten()) {
        assert!((x - y).abs() < 1e-12 + 1e-9 * y.abs());
    }
    let weight = STEEL.density * beam.members[0].section.area * 3.0 * 9.81;
    close(d.reaction[0][2], weight, 1e-9);
}

#[test]
fn a_mechanism_is_refused() {
    let mut beam = line([0.0; 3], [3.0, 0.0, 0.0], 2, Section::rectangle(0.1, 0.2));
    beam.supports = vec![(0, PINNED)];
    let case = Case { nodal: vec![(2, [0.0, 0.0, -1.0, 0.0, 0.0, 0.0])], ..Default::default() };
    assert!(solve(&beam, &case).unwrap_err().contains("mechanism"));
}

#[test]
fn solving_reports_its_stages() {
    let mut beam = line([0.0; 3], [3.0, 0.0, 0.0], 3, Section::rectangle(0.1, 0.2));
    beam.supports = vec![(0, FIXED)];
    let seen = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let into = seen.clone();
    let listening = sims_progress::listen_here(move |r| into.lock().unwrap().push(r.stage.clone()));
    solve(&beam, &Case { gravity: Some([0.0, 0.0, -9.81]), ..Default::default() }).unwrap();
    modes(&beam, 2).unwrap();
    drop(listening);
    let seen = seen.lock().unwrap();
    assert!(seen.iter().any(|s| s == "solving the frame"), "{seen:?}");
    assert!(seen.iter().any(|s| s == "finding the frame's modes"), "{seen:?}");
}

fn three_storey() -> (Model, impl Fn(usize, usize, usize) -> usize) {
    let (h, bx, by, storeys) = (3.2, 6.0, 5.0, 3);
    let concrete = Material::isotropic(30e9, 0.2, 0.0);
    let at = |i: usize, j: usize, k: usize| (k * 3 + j) * 3 + i;
    let mut model = Model::default();
    for k in 0..=storeys {
        for j in 0..3 {
            for i in 0..3 {
                model.nodes.push([i as f64 * bx, j as f64 * by, k as f64 * h]);
                if k == 0 {
                    model.supports.push((at(i, j, k), FIXED));
                } else {
                    model.masses.push((at(i, j, k), 20000.0));
                }
            }
        }
    }
    let bar =
        |a, b, section, up| Member { a, b, section, material: concrete, up, pinned: [false; 2] };
    for k in 0..storeys {
        for j in 0..3 {
            for i in 0..3 {
                model.members.push(bar(
                    at(i, j, k),
                    at(i, j, k + 1),
                    Section::rectangle(0.4, 0.4),
                    [1.0, 0.0, 0.0],
                ));
            }
        }
    }
    for k in 1..=storeys {
        for j in 0..3 {
            for i in 0..2 {
                model.members.push(bar(
                    at(i, j, k),
                    at(i + 1, j, k),
                    Section::rectangle(0.3, 0.5),
                    [0.0, 0.0, 1.0],
                ));
            }
        }
        for j in 0..2 {
            for i in 0..3 {
                model.members.push(bar(
                    at(i, j, k),
                    at(i, j + 1, k),
                    Section::rectangle(0.25, 0.45),
                    [0.0, 0.0, 1.0],
                ));
            }
        }
    }
    (model, at)
}

#[test]
fn a_three_storey_frame_matches_opensees() {
    let (model, at) = three_storey();
    let push = Case {
        nodal: (1..=3)
            .map(|k| (at(0, 0, k), [10000.0 * k as f64, 0.0, 0.0, 0.0, 0.0, 0.0]))
            .collect(),
        ..Default::default()
    };
    let r = solve(&model, &push).unwrap();
    close(r.displacement[at(0, 0, 3)][0], 0.0035598863289532947, 1e-6);
    close(r.displacement[at(2, 2, 3)][0], 0.0001314665213851835, 1e-6);
    close(r.displacement[at(1, 1, 3)][5], 0.00016878967530857794, 1e-6);
    close(r.reaction[at(0, 0, 0)][4], -29168.197563310197, 1e-6);
    let periods: Vec<f64> = modes(&model, 6).unwrap().iter().map(|m| 1.0 / m.frequency()).collect();
    let opensees = [
        0.670766630743724,
        0.6281759855060884,
        0.6194004766961124,
        0.4679225790556696,
        0.37203101867103155,
        0.3551936591838605,
    ];
    for (got, want) in periods.iter().zip(opensees) {
        close(*got, want, 1e-6);
    }
}
