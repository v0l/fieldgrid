use crate::bvh::Bvh;
use crate::grid::nearest;
use rayon::prelude::*;
use std::collections::HashSet;
use std::ops::RangeInclusive;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Filled {
    pub inside: Vec<usize>,
    pub thin: Vec<usize>,
}

fn spacing(line: &[f64]) -> f64 {
    line.windows(2).map(|w| w[1] - w[0]).filter(|d| *d > 0.0).fold(f64::INFINITY, f64::min)
}

fn span(line: &[f64], low: f64, high: f64) -> Option<RangeInclusive<usize>> {
    let last_sample = line.len().checked_sub(1)?;
    let first = line.partition_point(|s| *s <= low).saturating_sub(1).saturating_sub(2);
    let last = (line.partition_point(|s| *s < high) + 2).min(last_sample);
    (first <= last).then_some(first..=last)
}

pub fn fill(bvh: &Bvh, samples: [&[f64]; 3], index: impl Fn([usize; 3]) -> usize + Sync) -> Filled {
    let Some((low, high)) = bvh.bounds() else {
        return Filled::default();
    };
    let step = samples.map(|line| {
        let s = spacing(line);
        if s.is_finite() { s } else { (high[0] - low[0]).max(high[1] - low[1]).max(1.0) }
    });
    let mut spans = Vec::with_capacity(3);
    for a in 0..3 {
        let Some(s) = span(samples[a], low[a], high[a]) else {
            return Filled::default();
        };
        spans.push(s);
    }
    let mut inside = Vec::new();
    let mut thin = Vec::new();
    for axis in 0..3 {
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        let line = samples[axis];
        let n = line.len();
        let (us, vs) = (spans[u].clone(), spans[v].clone());
        let wide = us.end() - us.start() + 1;
        let tall = vs.end() - vs.start() + 1;
        let start = line[0].min(low[axis]) - 2.0 * step[axis];
        let mut direction = [0.0; 3];
        direction[axis] = 1.0;
        let reach = (line[0] - step[axis] * 0.5, line[n - 1] + step[axis] * 0.5);
        let columns: Vec<(Vec<usize>, Vec<usize>)> = (0..wide * tall)
            .into_par_iter()
            .map(|column| {
                let mut cell = [0usize; 3];
                cell[u] = us.start() + column % wide;
                cell[v] = vs.start() + column / wide;
                let mut from = [0.0; 3];
                from[u] = samples[u][cell[u]] + step[u] * 1.234_567e-6;
                from[v] = samples[v][cell[v]] + step[v] * 2.345_678e-6;
                from[axis] = start;
                let mut crossings = Vec::new();
                let mut last = 0.0;
                loop {
                    let t = bvh.ray(from, direction, usize::MAX, last + step[axis] * 1.0e-9);
                    if !t.is_finite() {
                        break;
                    }
                    crossings.push(start + t);
                    last = t;
                }
                let at = |k: usize| {
                    let mut c = cell;
                    c[axis] = k;
                    index(c)
                };
                let mut solid = Vec::new();
                let mut skins = Vec::new();
                for pair in crossings.chunks(2) {
                    let (a, b) = (pair[0], *pair.get(1).unwrap_or(&f64::INFINITY));
                    let first = line.partition_point(|s| *s < a);
                    let last = line.partition_point(|s| *s <= b).checked_sub(1);
                    match last {
                        Some(last) if first <= last => {
                            if axis == 0 {
                                solid.extend((first..=last).map(at));
                            }
                        }
                        _ if b.is_finite() => {
                            let middle = 0.5 * (a + b);
                            if middle >= reach.0 && middle < reach.1 {
                                skins.push(at(nearest(line, middle)));
                            }
                        }
                        _ => {}
                    }
                }
                (solid, skins)
            })
            .collect();
        for (solid, skins) in columns {
            inside.extend(solid);
            thin.extend(skins);
        }
    }
    let cored: HashSet<usize> = inside.iter().copied().collect();
    let mut thin: Vec<usize> = thin.into_iter().filter(|c| !cored.contains(c)).collect();
    thin.sort_unstable();
    thin.dedup();
    Filled { inside, thin }
}

pub type Triangle = [[f64; 3]; 3];

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

pub fn closed(triangles: &[Triangle]) -> bool {
    if triangles.is_empty() {
        return false;
    }
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in triangles.iter().flatten() {
        for a in 0..3 {
            lo[a] = lo[a].min(p[a]);
            hi[a] = hi[a].max(p[a]);
        }
    }
    let size = (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max).max(1e-30);
    let key =
        |p: [f64; 3]| p.map(|v| ((v - lo[0].min(lo[1]).min(lo[2])) / size * 1e9).round() as i64);
    let mut edges: std::collections::HashMap<([i64; 3], [i64; 3]), u32> = Default::default();
    for t in triangles {
        let k = t.map(key);
        if k[0] == k[1] || k[1] == k[2] || k[0] == k[2] {
            continue;
        }
        for (a, b) in [(k[0], k[1]), (k[1], k[2]), (k[2], k[0])] {
            *edges.entry(if a < b { (a, b) } else { (b, a) }).or_default() += 1;
        }
    }
    !edges.is_empty() && edges.values().all(|n| n % 2 == 0)
}

pub fn solid_edges(bvh: &Bvh, grid: &crate::grid::Grid) -> [Vec<bool>; 3] {
    let n = grid.dims();
    let total = n[0] * n[1] * n[2];
    std::array::from_fn(|c| {
        let centres = grid.centres(c);
        let samples: [&[f64]; 3] =
            std::array::from_fn(|a| if a == c { centres.as_slice() } else { grid.axis(a) });
        let filled = fill(bvh, samples, |p| (p[0] * n[1] + p[1]) * n[2] + p[2]);
        let mut mask = vec![false; total];
        for id in filled.inside.into_iter().chain(filled.thin) {
            mask[id] = true;
        }
        mask
    })
}

pub fn solid_cells(bvh: &Bvh, grid: &crate::grid::Grid) -> Vec<bool> {
    let n = grid.dims().map(|v| v - 1);
    let centres: [Vec<f64>; 3] = std::array::from_fn(|a| grid.centres(a));
    let samples: [&[f64]; 3] = std::array::from_fn(|a| centres[a].as_slice());
    let filled = fill(bvh, samples, |p| (p[0] * n[1] + p[1]) * n[2] + p[2]);
    let mut mask = vec![false; n[0] * n[1] * n[2]];
    for id in filled.inside.into_iter().chain(filled.thin) {
        mask[id] = true;
    }
    mask
}

fn local(line: &[f64], i: usize, centre: bool) -> f64 {
    let n = line.len();
    if centre {
        return line[(i + 1).min(n - 1)] - line[i.min(n - 2)];
    }
    let before = if i > 0 { line[i] - line[i - 1] } else { line[1] - line[0] };
    let after = if i + 1 < n { line[i + 1] - line[i] } else { line[n - 1] - line[n - 2] };
    0.5 * (before + after)
}

pub fn sheet_edges(triangles: &[Triangle], grid: &crate::grid::Grid) -> [Vec<bool>; 3] {
    let n = grid.dims();
    let total = n[0] * n[1] * n[2];
    let mut mask: [Vec<bool>; 3] = std::array::from_fn(|_| vec![false; total]);
    let lines = [grid.axis(0), grid.axis(1), grid.axis(2)];
    for t in triangles {
        let normal = cross(sub(t[1], t[0]), sub(t[2], t[0]));
        let area2 = dot(normal, normal).sqrt();
        if area2 <= 1e-30 {
            continue;
        }
        let nrm = normal.map(|v| v / area2);
        let lo: [f64; 3] =
            std::array::from_fn(|a| t.iter().map(|p| p[a]).fold(f64::INFINITY, f64::min));
        let hi: [f64; 3] =
            std::array::from_fn(|a| t.iter().map(|p| p[a]).fold(f64::NEG_INFINITY, f64::max));
        let range = |a: usize, centre: bool| {
            let line = lines[a];
            let cell = |x: f64| {
                let i = crate::grid::cell_of(line, x);
                line[i + 1] - line[i]
            };
            let reach = cell(lo[a]).max(cell(hi[a]));
            let first = line.partition_point(|v| *v < lo[a] - reach);
            let last = line.partition_point(|v| *v <= hi[a] + reach);
            let last = if centre { last.min(line.len() - 1) } else { last };
            first..last
        };
        let (e0, e1) = (sub(t[1], t[0]), sub(t[2], t[0]));
        let (d00, d01, d11) = (dot(e0, e0), dot(e0, e1), dot(e1, e1));
        let det = d00 * d11 - d01 * d01;
        for (c, edges) in mask.iter_mut().enumerate() {
            if nrm[c].abs() >= 0.75 {
                continue;
            }
            let ranges: [std::ops::Range<usize>; 3] = std::array::from_fn(|a| range(a, a == c));
            for i in ranges[0].clone() {
                for j in ranges[1].clone() {
                    for k in ranges[2].clone() {
                        let p = [i, j, k];
                        if p[c] + 1 >= n[c] {
                            continue;
                        }
                        let mid: [f64; 3] = std::array::from_fn(|a| {
                            if a == c {
                                0.5 * (lines[a][p[a]] + lines[a][p[a] + 1])
                            } else {
                                lines[a][p[a]]
                            }
                        });
                        let h: [f64; 3] = std::array::from_fn(|a| local(lines[a], p[a], a == c));
                        let tol = 0.5 * (0..3).map(|a| nrm[a].abs() * h[a]).sum::<f64>();
                        let d = dot(sub(mid, t[0]), nrm);
                        let eps = 1e-9 * tol;
                        if d <= eps - tol || d > tol + eps {
                            continue;
                        }
                        let q = sub(sub(mid, nrm.map(|v| v * d)), t[0]);
                        let (q0, q1) = (dot(q, e0), dot(q, e1));
                        let v = (d11 * q0 - d01 * q1) / det;
                        let w = (d00 * q1 - d01 * q0) / det;
                        let slack = 1e-9;
                        if v >= -slack && w >= -slack && v + w <= 1.0 + slack {
                            edges[(i * n[1] + j) * n[2] + k] = true;
                        }
                    }
                }
            }
        }
    }
    mask
}
