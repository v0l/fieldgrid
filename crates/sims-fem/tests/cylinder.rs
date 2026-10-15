use sims_fem::Quadratic;
use sims_fem::elastic::{Hold, Support, projectors};
use sims_tet::geom::*;
use sims_tet::{Point3, SurfaceMesh, Vector2, Vector3, fill};

const RADIUS: f64 = 2.5;
const LENGTH: f64 = 22.0;

#[derive(Clone, Copy)]
enum Shape {
    Wall,
    Bottom,
    Top,
}

impl Surface for Shape {
    fn at(&self, uv: Vector2) -> Option<Point3> {
        match self {
            Shape::Wall => Some(Point3::new(RADIUS * uv.x.cos(), RADIUS * uv.x.sin(), uv.y)),
            Shape::Bottom => Some(Point3::new(uv.x, uv.y, 0.0)),
            Shape::Top => Some(Point3::new(uv.x, uv.y, LENGTH)),
        }
    }

    fn nearest(&self, p: Point3, _: Option<Vector2>) -> Option<Vector2> {
        match self {
            Shape::Wall => Some(Vector2::new(p.y.atan2(p.x), p.z)),
            _ => Some(Vector2::new(p.x, p.y)),
        }
    }
}

fn cylinder(cell: f64) -> SurfaceMesh<Shape> {
    let rings = (std::f64::consts::TAU * RADIUS / cell).ceil() as usize;
    let rows = (LENGTH / cell).ceil() as usize;
    let mut mesh = SurfaceMesh {
        points: Vec::new(),
        triangles: Vec::new(),
        tags: Vec::new(),
        uvs: Vec::new(),
        faces: vec![Shape::Wall, Shape::Bottom, Shape::Top],
    };
    let mut index: std::collections::HashMap<u32, u32> = Default::default();
    let mut wall = |ring: usize, row: usize, mesh: &mut SurfaceMesh<Shape>| -> u32 {
        let key = (ring % rings) as u32 * 4096 + row as u32;
        *index.entry(key).or_insert_with(|| {
            let theta = std::f64::consts::TAU * (ring % rings) as f64 / rings as f64;
            let height = LENGTH * row as f64 / rows as f64;
            mesh.points.push(Point3::new(RADIUS * theta.cos(), RADIUS * theta.sin(), height));
            mesh.uvs.push(vec![(0, Vector2::new(theta, height))]);
            (mesh.points.len() - 1) as u32
        })
    };
    let mut triangles = Vec::new();
    for row in 0..rows {
        for ring in 0..rings {
            let q = [
                wall(ring, row, &mut mesh),
                wall(ring + 1, row, &mut mesh),
                wall(ring + 1, row + 1, &mut mesh),
                wall(ring, row + 1, &mut mesh),
            ];
            triangles.push(([q[0], q[1], q[2]], 0));
            triangles.push(([q[0], q[2], q[3]], 0));
        }
    }
    for (tag, level) in [(1usize, 0.0f64), (2, LENGTH)] {
        let mut flat: std::collections::HashMap<[i64; 2], u32> = Default::default();
        let mut point = |x: f64, y: f64, mesh: &mut SurfaceMesh<Shape>| -> u32 {
            let key = [(x * 1e6) as i64, (y * 1e6) as i64];
            *flat.entry(key).or_insert_with(|| {
                mesh.points.push(Point3::new(x, y, level));
                mesh.uvs.push(vec![(tag, Vector2::new(x, y))]);
                (mesh.points.len() - 1) as u32
            })
        };
        let spokes = rings;
        let ring_of = (RADIUS / cell).ceil() as usize;
        let mut at = |k: usize, mesh: &mut SurfaceMesh<Shape>| -> Vec<u32> {
            (0..=ring_of)
                .map(|r| {
                    let theta = std::f64::consts::TAU * (k % spokes) as f64 / spokes as f64;
                    let share = r as f64 / ring_of as f64;
                    point(RADIUS * share * theta.cos(), RADIUS * share * theta.sin(), mesh)
                })
                .collect()
        };
        for k in 0..spokes {
            let a = at(k, &mut mesh);
            let b = at(k + 1, &mut mesh);
            let (a, b) = if tag == 1 { (b, a) } else { (a, b) };
            for r in 0..a.len() - 1 {
                if a[r] == b[r] {
                    triangles.push(([a[r], b[r + 1], a[r + 1]], tag));
                    continue;
                }
                triangles.push(([a[r], b[r], b[r + 1]], tag));
                triangles.push(([a[r], b[r + 1], a[r + 1]], tag));
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
fn a_slide_on_a_curved_wall_holds_only_its_normal() {
    let tets = fill(cylinder(1.5), 1.5).expect("meshes");
    let mesh = Quadratic::from(&tets);
    let axis = Vector3::unit_z();
    for (tag, node, n) in &mesh.normals {
        if *tag != 0 {
            continue;
        }
        let lean = n.dot(axis).abs();
        assert!(lean < 0.02, "wall node {node} normal {n:?} leans {lean} on the axis");
    }
    assert!(!mesh.normals.is_empty(), "no surface normals recorded");
    let (blocks, owner) =
        projectors(&mesh, &[Support { label: "bore".into(), faces: vec![0], hold: Hold::Slide }]);
    let held = owner.iter().filter(|o| **o == Some(0)).count();
    assert!(held > 10, "the slide held {held} nodes");
    for (v, block) in blocks.iter().enumerate() {
        if owner[v] != Some(0) {
            continue;
        }
        let along = |d: [f64; 3]| {
            (0..3)
                .map(|k| (0..3).map(|j| block[k * 3 + j] * d[j]).sum::<f64>())
                .collect::<Vec<f64>>()
                .iter()
                .zip(d)
                .map(|(p, d)| p * d)
                .sum::<f64>()
        };
        let axial = along([0.0, 0.0, 1.0]);
        let p = mesh.nodes[v];
        let radial = [p.x / RADIUS, p.y / RADIUS, 0.0];
        let keeps = along(radial);
        assert!(
            (axial - 1.0).abs() < 0.02,
            "node {v} at {p:?} block {block:?} holds {:.3} of the axis; normals {:?}",
            1.0 - axial,
            mesh.normals
                .iter()
                .filter(|(t, n, _)| *t == 0 && *n as usize == v)
                .map(|(_, _, n)| n)
                .collect::<Vec<_>>()
        );
        assert!(
            keeps.abs() < 0.02,
            "node {v} at {p:?} block {block:?} frees {keeps} of the radial; normals {:?}",
            mesh.normals
                .iter()
                .filter(|(t, n, _)| *t == 0 && *n as usize == v)
                .map(|(_, _, n)| n)
                .collect::<Vec<_>>()
        );
    }
}
