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
