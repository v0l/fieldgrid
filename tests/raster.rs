use fieldgrid::Bvh;
use fieldgrid::grid::{Grading, centres, graded, uniform};
use fieldgrid::raster::fill;

fn cuboid(lo: [f64; 3], hi: [f64; 3]) -> Bvh {
    let p = |i: usize| [0, 1, 2].map(|a| if i >> a & 1 == 1 { hi[a] } else { lo[a] });
    let points: Vec<[f64; 3]> = (0..8).map(p).collect();
    let faces =
        [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]];
    let triangles: Vec<[usize; 3]> =
        faces.iter().flat_map(|f| [[f[0], f[1], f[2]], [f[0], f[2], f[3]]]).collect();
    Bvh::indexed(&points, &triangles)
}

fn sphere(radius: f64, rings: usize) -> Bvh {
    let segments = 2 * rings;
    let mut points = vec![[0.0, 0.0, radius], [0.0, 0.0, -radius]];
    for r in 1..rings {
        let theta = std::f64::consts::PI * r as f64 / rings as f64;
        for s in 0..segments {
            let phi = std::f64::consts::TAU * s as f64 / segments as f64;
            points.push([
                radius * theta.sin() * phi.cos(),
                radius * theta.sin() * phi.sin(),
                radius * theta.cos(),
            ]);
        }
    }
    let at = |r: usize, s: usize| 2 + (r - 1) * segments + s % segments;
    let mut triangles = Vec::new();
    for s in 0..segments {
        triangles.push([0, at(1, s), at(1, s + 1)]);
        triangles.push([1, at(rings - 1, s + 1), at(rings - 1, s)]);
        for r in 1..rings - 1 {
            triangles.push([at(r, s), at(r + 1, s), at(r + 1, s + 1)]);
            triangles.push([at(r, s), at(r + 1, s + 1), at(r, s + 1)]);
        }
    }
    Bvh::indexed(&points, &triangles)
}

fn xfast(n: [usize; 3]) -> impl Fn([usize; 3]) -> usize + Sync {
    move |[i, j, k]| i + (j + k * n[1]) * n[0]
}

fn unpack(n: [usize; 3], c: usize) -> [usize; 3] {
    [c % n[0], c / n[0] % n[1], c / (n[0] * n[1])]
}

#[test]
fn a_box_fills_exactly_the_samples_inside_it() {
    let line = uniform(0.0, 1.0, 20);
    let n = [line.len(); 3];
    let filled =
        fill(&cuboid([0.22, 0.31, 0.405], [0.71, 0.66, 0.9]), [&line, &line, &line], xfast(n));
    let want = |a: f64, b: f64| line.iter().filter(|s| **s > a && **s < b).count();
    assert_eq!(filled.inside.len(), want(0.22, 0.71) * want(0.31, 0.66) * want(0.405, 0.9));
    for c in &filled.inside {
        let [i, j, k] = unpack(n, *c);
        assert!(line[i] > 0.22 && line[i] < 0.71);
        assert!(line[j] > 0.31 && line[j] < 0.66);
        assert!(line[k] > 0.405 && line[k] < 0.9);
    }
    assert!(filled.thin.is_empty());
}

#[test]
fn a_sphere_fills_its_volume_on_a_graded_grid() {
    let line = graded(
        &[-1.0, 1.0],
        &Grading { features: &[(0.0, 0.01)], coarse: &|_| 0.05, ratio: 1.1, fill: 0.0 },
    );
    let cells = centres(&line);
    let widths: Vec<f64> = line.windows(2).map(|w| w[1] - w[0]).collect();
    let n = [cells.len(); 3];
    let filled = fill(&sphere(0.6, 48), [&cells, &cells, &cells], xfast(n));
    let volume: f64 = filled
        .inside
        .iter()
        .map(|c| unpack(n, *c).iter().map(|i| widths[*i]).product::<f64>())
        .sum();
    let exact = 4.0 / 3.0 * std::f64::consts::PI * 0.6f64.powi(3);
    assert!((volume / exact - 1.0).abs() < 0.02, "{volume} against {exact}");
}

#[test]
fn a_plate_thinner_than_a_cell_is_kept_as_a_skin() {
    let line = uniform(0.0, 1.0, 10);
    let n = [line.len(); 3];
    let filled =
        fill(&cuboid([0.15, 0.15, 0.52], [0.85, 0.85, 0.54]), [&line, &line, &line], xfast(n));
    assert!(filled.inside.is_empty());
    let k = line.iter().position(|z| (z - 0.5).abs() < 1e-12).unwrap();
    let inner = line.iter().filter(|s| **s > 0.15 && **s < 0.85).count();
    assert_eq!(filled.thin.len(), inner * inner);
    assert!(filled.thin.iter().all(|c| unpack(n, *c)[2] == k));
}

#[test]
fn nothing_to_fill_from_an_empty_mesh() {
    let line = uniform(0.0, 1.0, 4);
    let filled = fill(&Bvh::new(Vec::new()), [&line, &line, &line], |c| c[0]);
    assert!(filled.inside.is_empty() && filled.thin.is_empty());
}
