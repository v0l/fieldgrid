use fieldgrid::Bvh;
use fieldgrid::grid::{Grid, uniform};
use fieldgrid::raster::{Triangle, closed, sheet_edges, solid_cells, solid_edges};

fn quad(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> [Triangle; 2] {
    [[a, b, c], [a, c, d]]
}

fn plate(lo: [f64; 2], hi: [f64; 2], z: f64) -> Vec<Triangle> {
    quad([lo[0], lo[1], z], [hi[0], lo[1], z], [hi[0], hi[1], z], [lo[0], hi[1], z]).to_vec()
}

fn cuboid(lo: [f64; 3], hi: [f64; 3]) -> Vec<Triangle> {
    let p = |i: usize| [0, 1, 2].map(|a| if i >> a & 1 == 1 { hi[a] } else { lo[a] });
    let faces =
        [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]];
    faces.iter().flat_map(|f| quad(p(f[0]), p(f[1]), p(f[2]), p(f[3]))).collect()
}

fn grid() -> Grid {
    Grid { x: uniform(0.0, 1.0, 20), y: uniform(0.0, 1.0, 20), z: uniform(0.0, 1.0, 20), pml: 0 }
}

fn idx(n: [usize; 3], p: [usize; 3]) -> usize {
    (p[0] * n[1] + p[1]) * n[2] + p[2]
}

#[test]
fn a_box_is_closed_and_a_plate_is_not() {
    assert!(closed(&cuboid([0.1; 3], [0.4; 3])));
    assert!(!closed(&plate([0.1, 0.1], [0.4, 0.4], 0.5)));
    assert!(!closed(&[]));
}

#[test]
fn a_plate_on_a_grid_plane_marks_the_edges_it_covers() {
    let g = grid();
    let n = g.dims();
    let mask = sheet_edges(&plate([0.2, 0.3], [0.6, 0.5], 0.5), &g);
    for (c, edges) in mask.iter().enumerate() {
        for i in 0..n[0] {
            for j in 0..n[1] {
                for k in 0..n[2] {
                    let p = [i, j, k];
                    let x = |a: usize| g.axis(a)[p[a]];
                    let want = c != 2
                        && k == 10
                        && p[c] + 1 < n[c]
                        && (0..2).all(|a| {
                            let (lo, hi) = ([0.2, 0.3][a], [0.6, 0.5][a]);
                            if a == c {
                                x(a) >= lo - 1e-9 && g.axis(a)[p[a] + 1] <= hi + 1e-9
                            } else {
                                x(a) >= lo - 1e-9 && x(a) <= hi + 1e-9
                            }
                        });
                    assert_eq!(edges[idx(n, p)], want, "comp {c} at {p:?}");
                }
            }
        }
    }
}

#[test]
fn a_plate_between_planes_lands_on_the_nearest_one() {
    let g = grid();
    let n = g.dims();
    for (z, k) in [(0.512, 10), (0.538, 11), (0.525, 11), (0.5, 10)] {
        let mask = sheet_edges(&plate([0.2, 0.2], [0.6, 0.6], z), &g);
        let planes: std::collections::BTreeSet<usize> =
            (0..n[0] * n[1] * n[2]).filter(|id| mask[0][*id]).map(|id| id % n[2]).collect();
        assert_eq!(planes.into_iter().collect::<Vec<_>>(), vec![k], "z {z}");
    }
}

#[test]
fn a_solid_box_fills_the_edges_inside_it() {
    let g = grid();
    let n = g.dims();
    let mesh = cuboid([0.22, 0.31, 0.41], [0.68, 0.59, 0.77]);
    let bvh = Bvh::new(mesh.clone());
    let mask = solid_edges(&bvh, &g);
    let inside =
        |p: [f64; 3]| (0..3).all(|a| p[a] > [0.22, 0.31, 0.41][a] && p[a] < [0.68, 0.59, 0.77][a]);
    for (c, edges) in mask.iter().enumerate() {
        for i in 0..n[0] {
            for j in 0..n[1] {
                for k in 0..n[2] {
                    let p = [i, j, k];
                    if p[c] + 1 >= n[c] {
                        continue;
                    }
                    let mid: [f64; 3] = std::array::from_fn(|a| {
                        let l = g.axis(a);
                        if a == c { 0.5 * (l[p[a]] + l[p[a] + 1]) } else { l[p[a]] }
                    });
                    assert_eq!(edges[idx(n, p)], inside(mid), "comp {c} at {p:?}");
                }
            }
        }
    }
    let cells = solid_cells(&bvh, &g);
    let count = cells.iter().filter(|v| **v).count();
    let centre = |lo: f64, hi: f64| {
        (0..20)
            .filter(|i| {
                let x = (*i as f64 + 0.5) / 20.0;
                x > lo && x < hi
            })
            .count()
    };
    assert_eq!(count, centre(0.22, 0.68) * centre(0.31, 0.59) * centre(0.41, 0.77));
}

#[test]
fn a_slanted_plate_marks_a_staircase_that_stays_on_it() {
    let g = grid();
    let n = g.dims();
    let tri = quad([0.2, 0.2, 0.3], [0.7, 0.2, 0.3], [0.7, 0.6, 0.7], [0.2, 0.6, 0.7]);
    let mask = sheet_edges(&tri, &g);
    let mut count = 0;
    for (c, edges) in mask.iter().enumerate() {
        for id in (0..edges.len()).filter(|id| edges[*id]) {
            let p = [id / (n[1] * n[2]), id / n[2] % n[1], id % n[2]];
            let mid: [f64; 3] = std::array::from_fn(|a| {
                let l = g.axis(a);
                if a == c { 0.5 * (l[p[a]] + l[p[a] + 1]) } else { l[p[a]] }
            });
            let d = (mid[2] - 0.3) - (mid[1] - 0.2);
            assert!(d.abs() / 2f64.sqrt() <= 0.5 * 0.05 * 2f64.sqrt() + 1e-9, "{c} {p:?} {d}");
            count += 1;
        }
    }
    assert!(count > 200, "{count}");
}
