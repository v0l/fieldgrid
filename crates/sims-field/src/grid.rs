#[derive(Clone, Debug, Default, PartialEq)]
pub struct Grid {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
    pub pml: usize,
}

impl Grid {
    pub fn uniform(lo: [f64; 3], hi: [f64; 3], n: [usize; 3], pml: usize) -> Grid {
        let line = |a: usize| uniform(lo[a], hi[a], n[a]);
        Grid { x: line(0), y: line(1), z: line(2), pml }
    }

    pub fn axis(&self, a: usize) -> &[f64] {
        match a {
            0 => &self.x,
            1 => &self.y,
            _ => &self.z,
        }
    }

    pub fn dims(&self) -> [usize; 3] {
        [self.x.len(), self.y.len(), self.z.len()]
    }

    pub fn cells(&self) -> usize {
        (self.x.len() - 1) * (self.y.len() - 1) * (self.z.len() - 1)
    }

    pub fn nearest(&self, a: usize, v: f64) -> usize {
        nearest(self.axis(a), v)
    }

    pub fn cell_of(&self, a: usize, v: f64) -> usize {
        cell_of(self.axis(a), v)
    }

    pub fn centres(&self, a: usize) -> Vec<f64> {
        centres(self.axis(a))
    }
}

pub fn uniform(lo: f64, hi: f64, cells: usize) -> Vec<f64> {
    (0..=cells).map(|i| lo + (hi - lo) * i as f64 / cells as f64).collect()
}

pub fn nearest(line: &[f64], v: f64) -> usize {
    let i = line.partition_point(|x| *x < v);
    match (i.checked_sub(1), line.get(i)) {
        (Some(j), Some(x)) if (v - line[j]).abs() <= (x - v).abs() => j,
        (Some(j), None) => j,
        _ => i.min(line.len().saturating_sub(1)),
    }
}

pub fn cell_of(line: &[f64], v: f64) -> usize {
    line.partition_point(|x| *x <= v).saturating_sub(1).min(line.len().saturating_sub(2))
}

pub fn centres(line: &[f64]) -> Vec<f64> {
    line.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect()
}

pub struct Grading<'a> {
    pub features: &'a [(f64, f64)],
    pub coarse: &'a dyn Fn(f64) -> f64,
    pub ratio: f64,
    pub fill: f64,
}

pub fn graded(fixed: &[f64], grading: &Grading) -> Vec<f64> {
    let mut fixed = fixed.to_vec();
    fixed.sort_by(f64::total_cmp);
    fixed.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    let ratio = grading.ratio;
    let size = |x: f64| {
        grading
            .features
            .iter()
            .map(|(at, s)| s + (ratio - 1.0) * (x - at).abs())
            .fold((grading.coarse)(x), f64::min)
    };
    let Some(&first) = fixed.first() else {
        return Vec::new();
    };
    let mut out = vec![first];
    for w in fixed.windows(2) {
        let (a, b) = (w[0], w[1]);
        let steps = 256;
        let mut acc = vec![0.0; steps + 1];
        for k in 0..steps {
            let x = a + (b - a) * (k as f64 + 0.5) / steps as f64;
            acc[k + 1] = acc[k] + (b - a) / steps as f64 / size(x);
        }
        let fine = acc[steps].ceil().max(1.0) as usize;
        let n = if acc[steps] / fine as f64 >= grading.fill {
            fine
        } else {
            ((acc[steps] / grading.fill).floor() as usize).max(1)
        };
        for c in 1..n {
            let target = acc[steps] * c as f64 / n as f64;
            let k = acc.partition_point(|v| *v < target).clamp(1, steps);
            let t = (target - acc[k - 1]) / (acc[k] - acc[k - 1]).max(1e-30);
            out.push(a + (b - a) * ((k - 1) as f64 + t) / steps as f64);
        }
        out.push(b);
    }
    out
}

pub fn smooth(mut line: Vec<f64>, limit: f64) -> Vec<f64> {
    for _ in 0..200 {
        let d: Vec<f64> = line.windows(2).map(|w| w[1] - w[0]).collect();
        let bad = (0..d.len()).find(|&i| {
            (i > 0 && d[i] > limit * d[i - 1]) || (i + 1 < d.len() && d[i] > limit * d[i + 1])
        });
        let Some(i) = bad else {
            break;
        };
        let small = [i.checked_sub(1).map(|j| d[j]), d.get(i + 1).copied()]
            .into_iter()
            .flatten()
            .fold(f64::INFINITY, f64::min);
        let parts = (d[i] / (limit * small)).ceil().max(2.0) as usize;
        let (a, b) = (line[i], line[i + 1]);
        let extra: Vec<f64> = (1..parts).map(|k| a + (b - a) * k as f64 / parts as f64).collect();
        line.splice(i + 1..i + 1, extra);
    }
    line
}

pub fn pad(mut line: Vec<f64>, cells: usize) -> Vec<f64> {
    if line.len() < 2 {
        return line;
    }
    let first = line[1] - line[0];
    let last = line[line.len() - 1] - line[line.len() - 2];
    let mut pre: Vec<f64> = (1..=cells).rev().map(|i| line[0] - first * i as f64).collect();
    let post: Vec<f64> = (1..=cells).map(|i| line[line.len() - 1] + last * i as f64).collect();
    pre.append(&mut line);
    pre.extend(post);
    pre
}
