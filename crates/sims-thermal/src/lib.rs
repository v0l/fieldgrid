use faer::linalg::solvers::Solve;

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
        (p[2] * self.y.len() + p[1]) * self.x.len() + p[0]
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
    edges: Vec<(usize, usize, f64)>,
    films: Vec<(usize, usize, f64)>,
    fixed: Vec<Vec<usize>>,
    index: Vec<usize>,
    order: Vec<usize>,
    factor: Option<faer::sparse::linalg::solvers::Llt<usize, f64>>,
    environments: usize,
}

const NONE: usize = usize::MAX;

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
        let _stage = sims_progress::stage("factoring the conduction grid", None, "");
        let g = &self.grid;
        let n = g.nodes();
        if n.iter().any(|c| *c < 2) {
            return Err("the grid needs at least two lines along each axis".into());
        }
        let total = n[0] * n[1] * n[2];
        let mut edges: std::collections::HashMap<(usize, usize), f64> = Default::default();
        let cells = [n[0] - 1, n[1] - 1, n[2] - 1];
        let cell = |c: [usize; 3]| (c[2] * cells[1] + c[1]) * cells[0] + c[0];
        let mut lambdas = vec![0.0; cells[0] * cells[1] * cells[2]];
        for k in 0..n[2] - 1 {
            for j in 0..n[1] - 1 {
                for i in 0..n[0] - 1 {
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
                        let g4 = lambda * d[u] * d[v] / (4.0 * d[a]);
                        for su in 0..2 {
                            for sv in 0..2 {
                                let mut p = [i, j, k];
                                p[u] += su;
                                p[v] += sv;
                                let mut q = p;
                                q[a] += 1;
                                *edges.entry((g.index(p), g.index(q))).or_default() += g4;
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
        let mut connected = vec![false; total];
        for (a, b) in edges.keys() {
            connected[*a] = true;
            connected[*b] = true;
        }
        let mut index = vec![NONE; total];
        let mut order = Vec::new();
        for i in 0..total {
            if connected[i] && fixed_by[i].is_empty() {
                index[i] = order.len();
                order.push(i);
            }
        }
        let mut lower: std::collections::HashMap<(usize, usize), f64> = Default::default();
        for ((a, b), w) in &edges {
            for (p, q) in [(*a, *b), (*b, *a)] {
                if index[p] != NONE {
                    *lower.entry((index[p], index[p])).or_default() += w;
                    if index[q] != NONE && index[q] < index[p] {
                        *lower.entry((index[p], index[q])).or_default() -= w;
                    }
                }
            }
        }
        for ((node, _), w) in &films {
            if index[*node] != NONE {
                *lower.entry((index[*node], index[*node])).or_default() += w;
            }
        }
        let m = order.len();
        let factor = if m > 0 {
            let triplets: Vec<faer::sparse::Triplet<usize, usize, f64>> =
                lower.into_iter().map(|((r, c), v)| faer::sparse::Triplet::new(r, c, v)).collect();
            let a =
                faer::sparse::SparseColMat::<usize, f64>::try_new_from_triplets(m, m, &triplets)
                    .map_err(|e| format!("{e:?}"))?;
            Some(a.sp_cholesky(faer::Side::Lower).map_err(|_| {
                "part of the construction touches no boundary with a temperature".to_string()
            })?)
        } else {
            None
        };
        Ok(Prepared {
            grid: g.clone(),
            edges: edges.into_iter().map(|((a, b), w)| (a, b, w)).collect(),
            films: films.into_iter().map(|((n, e), w)| (n, e, w)).collect(),
            fixed: fixed_by,
            index,
            order,
            factor,
            environments,
        })
    }
}

impl Prepared {
    pub fn solve(&self, temperatures: &[f64]) -> Result<Solution, String> {
        if temperatures.len() < self.environments {
            return Err(format!(
                "the boundaries name {} environments but {} temperatures were given",
                self.environments,
                temperatures.len()
            ));
        }
        let _stage = sims_progress::stage("solving the conduction grid", None, "");
        let total = self.index.len();
        let mut t = vec![f64::NAN; total];
        for (i, envs) in self.fixed.iter().enumerate() {
            if !envs.is_empty() {
                t[i] = envs.iter().map(|e| temperatures[*e]).sum::<f64>() / envs.len() as f64;
            }
        }
        let m = self.order.len();
        if let Some(factor) = &self.factor {
            let mut b = faer::Mat::<f64>::zeros(m, 1);
            for (node, env, w) in &self.films {
                if self.index[*node] != NONE {
                    b[(self.index[*node], 0)] += w * temperatures[*env];
                }
            }
            for (p, q, w) in &self.edges {
                for (a, c) in [(*p, *q), (*q, *p)] {
                    if self.index[a] != NONE && !self.fixed[c].is_empty() {
                        b[(self.index[a], 0)] += w * t[c];
                    }
                }
            }
            let x = factor.solve(&b);
            for (k, node) in self.order.iter().enumerate() {
                t[*node] = x[(k, 0)];
            }
        }
        let mut flows = vec![0.0; self.environments.max(temperatures.len())];
        for (node, env, w) in &self.films {
            if !t[*node].is_nan() {
                flows[*env] += w * (temperatures[*env] - t[*node]);
            }
        }
        let mut out = vec![0.0; total];
        for (p, q, w) in &self.edges {
            if t[*p].is_nan() || t[*q].is_nan() {
                continue;
            }
            out[*p] += w * (t[*p] - t[*q]);
            out[*q] += w * (t[*q] - t[*p]);
        }
        for (i, envs) in self.fixed.iter().enumerate() {
            for e in envs {
                flows[*e] += out[i] / envs.len() as f64;
            }
        }
        Ok(Solution { temperature: t, flows })
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
