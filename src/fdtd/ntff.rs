use super::{C0, ETA0, Sim, idx};
use std::f64::consts::PI;

#[derive(Clone, Debug, Default)]
pub struct Patches {
    pub pos: Vec<[f64; 3]>,
    pub area: Vec<f64>,
    pub gpu: Vec<u32>,
}

pub fn patches(sim: &Sim, inset: usize) -> Patches {
    let n = sim.dims();
    let p = sim.pml + inset;
    let lo = [p; 3];
    let hi = [n[0] - 1 - p, n[1] - 1 - p, n[2] - 1 - p];
    let g = &sim.grid;
    let mut out = Patches { pos: Vec::new(), area: Vec::new(), gpu: Vec::new() };
    for axis in 0..3 {
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        for (side, plane) in [(-1i8, lo[axis]), (1i8, hi[axis])] {
            for a in lo[u]..hi[u] {
                for b in lo[v]..hi[v] {
                    let mut pos = [0.0; 3];
                    pos[axis] = g.axis(axis)[plane];
                    pos[u] = 0.5 * (g.axis(u)[a] + g.axis(u)[a + 1]);
                    pos[v] = 0.5 * (g.axis(v)[b] + g.axis(v)[b + 1]);
                    out.pos.push(pos);
                    out.area.push(sim.ax[u].d[a] * sim.ax[v].d[b]);
                    let mut q = [0usize; 3];
                    q[axis] = plane;
                    q[u] = a;
                    q[v] = b;
                    let at = |moves: &[(usize, isize)]| {
                        let mut r = [q[0] as isize, q[1] as isize, q[2] as isize];
                        for &(ax, d) in moves {
                            r[ax] += d;
                        }
                        idx(n, r[0] as usize, r[1] as usize, r[2] as usize) as u32
                    };
                    out.gpu.push(axis as u32 | if side > 0 { 4 } else { 0 });
                    out.gpu.extend([at(&[]), at(&[(v, 1)]), at(&[]), at(&[(u, 1)])]);
                    out.gpu.extend([
                        at(&[]),
                        at(&[(axis, -1)]),
                        at(&[(u, 1)]),
                        at(&[(axis, -1), (u, 1)]),
                    ]);
                    out.gpu.extend([
                        at(&[]),
                        at(&[(axis, -1)]),
                        at(&[(v, 1)]),
                        at(&[(axis, -1), (v, 1)]),
                    ]);
                }
            }
        }
    }
    out
}

pub type C = (f64, f64);

fn mul(a: C, b: C) -> C {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}

#[derive(Clone, Debug, Default)]
pub struct Far {
    pub pos: Vec<[f64; 3]>,
    pub area: Vec<f64>,
    pub j: Vec<[C; 3]>,
    pub m: Vec<[C; 3]>,
    pub freq: f64,
}

impl Far {
    pub fn load(patches: &Patches, acc: &[f32], nf: usize, k: usize, freq: f64, scale: C) -> Far {
        let count = patches.pos.len();
        let mut j = Vec::with_capacity(count);
        let mut m = Vec::with_capacity(count);
        for q in 0..count {
            let o = (q * nf + k) * 12;
            let c = |x: usize| mul((acc[o + x] as f64, acc[o + x + 1] as f64), scale);
            j.push([c(0), c(2), c(4)]);
            m.push([c(6), c(8), c(10)]);
        }
        Far { pos: patches.pos.clone(), area: patches.area.clone(), j, m, freq }
    }

    pub fn field(&self, dir: [f64; 3]) -> [C; 3] {
        let k = 2.0 * PI * self.freq / C0;
        let mut nn = [(0.0, 0.0); 3];
        let mut ll = [(0.0, 0.0); 3];
        for (i, p) in self.pos.iter().enumerate() {
            let arg = k * (dir[0] * p[0] + dir[1] * p[1] + dir[2] * p[2]);
            let ph = (self.area[i] * arg.cos(), self.area[i] * arg.sin());
            for c in 0..3 {
                let a = mul(self.j[i][c], ph);
                let b = mul(self.m[i][c], ph);
                nn[c] = (nn[c].0 + a.0, nn[c].1 + a.1);
                ll[c] = (ll[c].0 + b.0, ll[c].1 + b.1);
            }
        }
        let rn = (0..3).fold((0.0, 0.0), |s, c| (s.0 + nn[c].0 * dir[c], s.1 + nn[c].1 * dir[c]));
        let perp: [C; 3] =
            std::array::from_fn(|c| (nn[c].0 - rn.0 * dir[c], nn[c].1 - rn.1 * dir[c]));
        let rl = [
            (dir[1] * ll[2].0 - dir[2] * ll[1].0, dir[1] * ll[2].1 - dir[2] * ll[1].1),
            (dir[2] * ll[0].0 - dir[0] * ll[2].0, dir[2] * ll[0].1 - dir[0] * ll[2].1),
            (dir[0] * ll[1].0 - dir[1] * ll[0].0, dir[0] * ll[1].1 - dir[1] * ll[0].1),
        ];
        let scale = k / (4.0 * PI);
        std::array::from_fn(|c| {
            ((rl[c].0 - perp[c].0 * ETA0) * scale, (rl[c].1 - perp[c].1 * ETA0) * scale)
        })
    }

    pub fn intensity(&self, dir: [f64; 3]) -> f64 {
        self.field(dir).iter().map(|c| c.0 * c.0 + c.1 * c.1).sum::<f64>() / (2.0 * ETA0)
    }

    pub fn radiated(&self, steps: usize) -> (f64, f64) {
        let mut total = 0.0;
        let mut peak: f64 = 0.0;
        for a in 0..steps {
            let th = PI * (a as f64 + 0.5) / steps as f64;
            for b in 0..2 * steps {
                let ph = PI * b as f64 / steps as f64;
                let dir = [th.sin() * ph.cos(), th.sin() * ph.sin(), th.cos()];
                let u = self.intensity(dir);
                peak = peak.max(u);
                total += u * th.sin() * (PI / steps as f64) * (PI / steps as f64);
            }
        }
        (total, peak)
    }
}
