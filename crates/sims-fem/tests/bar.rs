use sims_fem::Quadratic;
use sims_fem::elastic::{Hold, Load, Problem, Push, Support, solve};
use sims_tet::{HashMap, Point3, SurfaceMesh, Vector3, fill};

const LENGTH: f64 = 100.0;
const SIDE: f64 = 10.0;

fn bar(cell: f64) -> SurfaceMesh<()> {
    let size = [LENGTH, SIDE, SIDE];
    let mut mesh = SurfaceMesh {
        points: Vec::new(),
        triangles: Vec::new(),
        tags: Vec::new(),
        uvs: Vec::new(),
        faces: vec![(); 6],
    };
    let mut index: HashMap<[i64; 3], u32> = HashMap::default();
    let mut point = |p: [f64; 3]| -> u32 {
        let key = p.map(|v| (v * 1e6).round() as i64);
        *index.entry(key).or_insert_with(|| {
            mesh.points.push(Point3::new(p[0], p[1], p[2]));
            mesh.uvs.push(Vec::new());
            (mesh.points.len() - 1) as u32
        })
    };
    let mut triangles = Vec::new();
    for axis in 0..3 {
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        let (nu, nv) = ((size[u] / cell).ceil() as usize, (size[v] / cell).ceil() as usize);
        for high in [false, true] {
            let tag = axis * 2 + high as usize;
            let at = |i: usize, j: usize| {
                let mut p = [0.0; 3];
                p[axis] = if high { size[axis] } else { 0.0 };
                p[u] = size[u] * i as f64 / nu as f64;
                p[v] = size[v] * j as f64 / nv as f64;
                p
            };
            for i in 0..nu {
                for j in 0..nv {
                    let q =
                        [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)].map(&mut point);
                    let (a, b) = if high {
                        ([q[0], q[1], q[2]], [q[0], q[2], q[3]])
                    } else {
                        ([q[0], q[2], q[1]], [q[0], q[3], q[2]])
                    };
                    triangles.push((a, tag));
                    triangles.push((b, tag));
                }
            }
        }
    }
    for (t, tag) in triangles {
        mesh.triangles.push(t);
        mesh.tags.push(tag);
    }
    mesh
}

#[test]
fn a_bar_in_tension_stretches_by_fl_over_ae() {
    let tets = fill(bar(2.5), 2.5).expect("meshes");
    let quality = tets.quality();
    assert!((quality.volume / (LENGTH * SIDE * SIDE) - 1.0).abs() < 1e-9, "{quality:?}");
    let mesh = Quadratic::from(&tets);
    let (young, force) = (200e3, 1000.0);
    let answer = solve(&Problem {
        mesh: &mesh,
        materials: &[(young, 0.0)],
        supports: &[Support { label: "end".into(), faces: vec![0], hold: Hold::Fixed }],
        loads: &[Load {
            label: "pull".into(),
            faces: vec![1],
            push: Push::Force(Vector3::new(force, 0.0, 0.0)),
        }],
        body_force: None,
        expansion: None,
        wind: None,
    })
    .expect("solves");
    let tip: Vec<f64> = mesh
        .nodes
        .iter()
        .zip(&answer.displacement)
        .filter(|(p, _)| (p.x - LENGTH).abs() < 1e-9)
        .map(|(_, d)| d.x)
        .collect();
    let mean = tip.iter().sum::<f64>() / tip.len() as f64;
    let exact = force * LENGTH / (young * SIDE * SIDE);
    assert!((mean / exact - 1.0).abs() < 0.01, "{mean} against {exact}");
    let reaction = answer.reactions.iter().find(|r| r.0 == "end").expect("a reaction").1;
    assert!((reaction.x + force).abs() < force * 1e-6, "{reaction:?}");
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
fn meshing_and_solving_report_progress() {
    let stages = stages_of(|| {
        let tets = fill(bar(5.0), 5.0).expect("meshes");
        let mesh = Quadratic::from(&tets);
        solve(&Problem {
            mesh: &mesh,
            materials: &[(200e3, 0.3)],
            supports: &[Support { label: "end".into(), faces: vec![0], hold: Hold::Fixed }],
            loads: &[Load {
                label: "pull".into(),
                faces: vec![1],
                push: Push::Force(Vector3::new(1000.0, 0.0, 0.0)),
            }],
            body_force: None,
            expansion: None,
            wind: None,
        })
        .expect("solves");
    });
    assert!(stages.iter().any(|s| s == "filling the volume"), "{stages:?}");
    assert!(stages.len() > 1, "{stages:?}");
}

fn bonded_nodes(gap: f64, apart: &[(usize, usize)]) -> usize {
    let left = Quadratic::from(&fill(bar(2.5), 2.5).expect("meshes"));
    let mut right = Quadratic::from(&fill(bar(2.5), 2.5).expect("meshes"));
    right.nodes.iter_mut().for_each(|p| p.x += LENGTH + gap);
    let parts = vec![("left".to_string(), left), ("right".to_string(), right)];
    let merged = sims_fem::merge(parts, &[0.25, 0.25], &[], &[], apart).expect("merges");
    merged.bonds.iter().map(|b| b.2).sum()
}

#[test]
fn bars_bond_across_a_small_gap_unless_kept_apart() {
    assert!(bonded_nodes(0.0, &[]) > 50);
    assert!(bonded_nodes(0.2, &[]) > 50);
    assert_eq!(bonded_nodes(0.2, &[(0, 1)]), 0);
    assert_eq!(bonded_nodes(0.4, &[]), 0);
}
