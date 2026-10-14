pub mod nodal;

use sims_field::Gpu;

#[derive(Clone, Debug, PartialEq)]
pub struct Grid {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
}

impl Grid {
    pub fn plane(x: Vec<f64>, y: Vec<f64>) -> Grid {
        Grid { x, y, z: vec![0.0, 1.0] }
    }

    pub fn nodes(&self) -> [usize; 3] {
        [self.x.len(), self.y.len(), self.z.len()]
    }

    fn axis(&self, a: usize) -> &[f64] {
        match a {
            0 => &self.x,
            1 => &self.y,
            _ => &self.z,
        }
    }

    fn index(&self, p: [usize; 3]) -> usize {
        (p[0] * self.y.len() + p[1]) * self.z.len() + p[2]
    }
}

pub fn lines(stops: &[f64], fine: f64, grow: f64) -> Vec<f64> {
    let mut stops: Vec<f64> = stops.to_vec();
    stops.sort_by(f64::total_cmp);
    stops.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    let mut out = vec![stops[0]];
    for w in stops.windows(2) {
        let (a, b) = (w[0], w[1]);
        let mut left = vec![];
        let mut right = vec![];
        let (mut lo, mut hi, mut step) = (a, b, fine);
        while hi - lo > 3.0 * step * (1.0 - 1e-9) {
            lo += step;
            hi -= step;
            left.push(lo);
            right.push(hi);
            step *= grow;
        }
        let gap = hi - lo;
        if gap > 1.5 * step {
            left.push(lo + gap / 2.0);
        }
        out.extend(left);
        out.extend(right.into_iter().rev());
        out.push(b);
    }
    out
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Block {
    pub low: [f64; 3],
    pub high: [f64; 3],
    pub conductivity: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Boundary {
    pub low: [f64; 3],
    pub high: [f64; 3],
    pub environment: usize,
    pub resistance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Model {
    pub grid: Grid,
    pub blocks: Vec<Block>,
    pub boundaries: Vec<Boundary>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Solution {
    pub temperature: Vec<f64>,
    pub flows: Vec<f64>,
}

pub struct Prepared {
    grid: Grid,
    problem: nodal::Problem,
    films: Vec<(usize, usize, f64)>,
    fixed: Vec<(usize, Vec<usize>)>,
    connected: Vec<bool>,
    environments: usize,
}

fn inside(p: f64, lo: f64, hi: f64) -> bool {
    p >= lo.min(hi) - 1e-12 && p <= lo.max(hi) + 1e-12
}

impl Model {
    fn conductivity(&self, centre: [f64; 3]) -> f64 {
        self.blocks
            .iter()
            .rev()
            .find(|b| (0..3).all(|a| inside(centre[a], b.low[a], b.high[a])))
            .map_or(0.0, |b| b.conductivity)
    }

    pub fn prepare(&self) -> Result<Prepared, String> {
        let g = &self.grid;
        let n = g.nodes();
        if n.iter().any(|c| *c < 2) {
            return Err("the grid needs at least two lines along each axis".into());
        }
        let total = n[0] * n[1] * n[2];
        let cells = [n[0] - 1, n[1] - 1, n[2] - 1];
        let cell = |c: [usize; 3]| (c[2] * cells[1] + c[1]) * cells[0] + c[0];
        let mut lambdas = vec![0.0; cells[0] * cells[1] * cells[2]];
        let mut along = [vec![0.0f64; total], vec![0.0f64; total], vec![0.0f64; total]];
        let mut connected = vec![false; total];
        for k in 0..cells[2] {
            for j in 0..cells[1] {
                for i in 0..cells[0] {
                    let d = [g.x[i + 1] - g.x[i], g.y[j + 1] - g.y[j], g.z[k + 1] - g.z[k]];
                    let centre = [
                        (g.x[i] + g.x[i + 1]) / 2.0,
                        (g.y[j] + g.y[j + 1]) / 2.0,
                        (g.z[k] + g.z[k + 1]) / 2.0,
                    ];
                    let lambda = self.conductivity(centre);
                    lambdas[cell([i, j, k])] = lambda;
                    if lambda <= 0.0 {
                        continue;
                    }
                    for a in 0..3 {
                        let (u, v) = ((a + 1) % 3, (a + 2) % 3);
                        let quarter = lambda * d[u] * d[v] / (4.0 * d[a]);
                        for su in 0..2 {
                            for sv in 0..2 {
                                let mut p = [i, j, k];
                                p[u] += su;
                                p[v] += sv;
                                along[a][g.index(p)] += quarter;
                                connected[g.index(p)] = true;
                                p[a] += 1;
                                connected[g.index(p)] = true;
                            }
                        }
                    }
                }
            }
        }
        let environments = self.boundaries.iter().map(|b| b.environment + 1).max().unwrap_or(0);
        let mut films: std::collections::HashMap<(usize, usize), f64> = Default::default();
        let mut fixed_by: Vec<Vec<usize>> = vec![Vec::new(); total];
        for k in 0..cells[2] {
            for j in 0..cells[1] {
                for i in 0..cells[0] {
                    let c = [i, j, k];
                    if lambdas[cell(c)] <= 0.0 {
                        continue;
                    }
                    for axis in 0..3 {
                        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
                        for high in [false, true] {
                            let open = if high { c[axis] + 1 >= cells[axis] } else { c[axis] == 0 };
                            let exposed = open || {
                                let mut nb = c;
                                if high {
                                    nb[axis] += 1
                                } else {
                                    nb[axis] -= 1
                                }
                                lambdas[cell(nb)] <= 0.0
                            };
                            if !exposed {
                                continue;
                            }
                            let plane = c[axis] + high as usize;
                            let mut centre = [0.0; 3];
                            for a in 0..3 {
                                let ax = g.axis(a);
                                centre[a] = if a == axis {
                                    ax[plane]
                                } else {
                                    (ax[c[a]] + ax[c[a] + 1]) / 2.0
                                };
                            }
                            let Some(b) =
                                self.boundaries.iter().rev().find(|b| {
                                    (0..3).all(|a| inside(centre[a], b.low[a], b.high[a]))
                                })
                            else {
                                continue;
                            };
                            let (au, av) = (g.axis(u), g.axis(v));
                            let area = (au[c[u] + 1] - au[c[u]]) * (av[c[v] + 1] - av[c[v]]);
                            for su in 0..2 {
                                for sv in 0..2 {
                                    let mut p = c;
                                    p[axis] = plane;
                                    p[u] += su;
                                    p[v] += sv;
                                    let node = g.index(p);
                                    if b.resistance > 0.0 {
                                        *films.entry((node, b.environment)).or_default() +=
                                            area / (4.0 * b.resistance);
                                    } else if !fixed_by[node].contains(&b.environment) {
                                        fixed_by[node].push(b.environment);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        let mut problem = nodal::Problem::new(n);
        let [gx, gy, gz] = along;
        problem.gx = gx.into_iter().map(|v| v as f32).collect();
        problem.gy = gy.into_iter().map(|v| v as f32).collect();
        problem.gz = gz.into_iter().map(|v| v as f32).collect();
        Ok(Prepared {
            grid: g.clone(),
            problem,
            films: films.into_iter().map(|((node, e), w)| (node, e, w)).collect(),
            fixed: fixed_by.into_iter().enumerate().filter(|(_, e)| !e.is_empty()).collect(),
            connected,
            environments,
        })
    }
}

impl Prepared {
    fn problem(&self, temperatures: &[f64]) -> Result<nodal::Problem, String> {
        if temperatures.len() < self.environments {
            return Err(format!(
                "the boundaries name {} environments but {} temperatures were given",
                self.environments,
                temperatures.len()
            ));
        }
        let mut p = nodal::Problem::new(self.problem.n);
        p.gx.clone_from(&self.problem.gx);
        p.gy.clone_from(&self.problem.gy);
        p.gz.clone_from(&self.problem.gz);
        let mut g0 = vec![0.0f64; p.len()];
        let mut source = vec![0.0f64; p.len()];
        for (node, env, w) in &self.films {
            g0[*node] += w;
            source[*node] += w * temperatures[*env];
        }
        p.g0 = g0.into_iter().map(|v| v as f32).collect();
        p.source = source.into_iter().map(|v| v as f32).collect();
        for (node, envs) in &self.fixed {
            let t = envs.iter().map(|e| temperatures[*e]).sum::<f64>() / envs.len() as f64;
            p.fixed[*node] = Some(t as f32);
        }
        Ok(p)
    }

    fn finish(
        &self,
        temperatures: &[f64],
        problem: &nodal::Problem,
        s: nodal::Solution,
    ) -> Solution {
        let t: Vec<f64> = (0..problem.len())
            .map(|i| match problem.fixed[i] {
                Some(v) => v as f64,
                None if self.connected[i] => s.delta[i],
                None => f64::NAN,
            })
            .collect();
        let mut flows = vec![0.0; self.environments.max(temperatures.len())];
        for (node, env, w) in &self.films {
            flows[*env] += w * (temperatures[*env] - t[*node]);
        }
        let [nx, ny, nz] = problem.n;
        let stride = [ny * nz, nz, 1];
        let dims = [nx, ny, nz];
        for (node, envs) in &self.fixed {
            let p = [node / stride[0], node / stride[1] % ny, node % nz];
            let mut out = 0.0;
            for (a, g) in [&problem.gx, &problem.gy, &problem.gz].into_iter().enumerate() {
                if p[a] + 1 < dims[a] {
                    out += g[*node] as f64 * (t[*node] - t[node + stride[a]]);
                }
                if p[a] > 0 {
                    out += g[node - stride[a]] as f64 * (t[*node] - t[node - stride[a]]);
                }
            }
            for e in envs {
                flows[*e] += out / envs.len() as f64;
            }
        }
        Solution { temperature: t, flows }
    }

    pub fn solve(&self, temperatures: &[f64]) -> Result<Solution, String> {
        let problem = self.problem(temperatures)?;
        let s = problem.solve_direct().map_err(|_| {
            "part of the construction touches no boundary with a temperature".to_string()
        })?;
        Ok(self.finish(temperatures, &problem, s))
    }

    pub fn solve_iterative(
        &self,
        temperatures: &[f64],
        gpu: Option<&Gpu>,
        tolerance: f64,
    ) -> Result<Solution, String> {
        let problem = self.problem(temperatures)?;
        let s = problem.solve(gpu, tolerance, 100 * problem.len().max(1000));
        if s.residual.is_nan() || s.residual > tolerance {
            return Err(format!("the solve stopped at a residual of {:.1e}", s.residual));
        }
        Ok(self.finish(temperatures, &problem, s))
    }

    pub fn coupling(&self) -> Result<Vec<Vec<f64>>, String> {
        let n = self.environments;
        (0..n)
            .map(|j| {
                let mut unit = vec![0.0; n];
                unit[j] = 1.0;
                let s = self.solve(&unit)?;
                Ok((0..n).map(|i| if i == j { s.flows[i] } else { -s.flows[i] }).collect())
            })
            .collect()
    }

    pub fn temperature(&self, solution: &Solution, p: [f64; 3]) -> Option<f64> {
        let g = &self.grid;
        let mut cell = [0usize; 3];
        let mut frac = [0.0; 3];
        for a in 0..3 {
            let ax = g.axis(a);
            if p[a] < ax[0] - 1e-12 || p[a] > ax[ax.len() - 1] + 1e-12 {
                return None;
            }
            let c = ax.partition_point(|v| *v <= p[a]).clamp(1, ax.len() - 1) - 1;
            cell[a] = c;
            frac[a] = ((p[a] - ax[c]) / (ax[c + 1] - ax[c])).clamp(0.0, 1.0);
        }
        let mut sum = 0.0;
        let mut weight = 0.0;
        for corner in 0..8 {
            let mut q = cell;
            let mut w = 1.0;
            for a in 0..3 {
                let up = corner >> a & 1 == 1;
                if up {
                    q[a] += 1;
                }
                w *= if up { frac[a] } else { 1.0 - frac[a] };
            }
            let v = solution.temperature[g.index(q)];
            if w > 0.0 && !v.is_nan() {
                sum += w * v;
                weight += w;
            }
        }
        (weight > 0.0).then(|| sum / weight)
    }
}
