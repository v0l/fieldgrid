#![allow(clippy::needless_range_loop)]

mod element;

use element::{Condensed, M12, V12};
use faer::linalg::solvers::Solve;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Section {
    pub area: f64,
    pub iy: f64,
    pub iz: f64,
    pub j: f64,
}

impl Section {
    pub fn rectangle(width: f64, depth: f64) -> Section {
        let (a, b) = (width.max(depth), width.min(depth));
        let j = a * b.powi(3) * (1.0 / 3.0 - 0.21 * b / a * (1.0 - b.powi(4) / (12.0 * a.powi(4))));
        Section {
            area: width * depth,
            iy: depth * width.powi(3) / 12.0,
            iz: width * depth.powi(3) / 12.0,
            j,
        }
    }

    pub fn tube(outer: f64, wall: f64) -> Section {
        let (ro, ri) = (outer / 2.0, outer / 2.0 - wall);
        let i = std::f64::consts::PI / 4.0 * (ro.powi(4) - ri.powi(4));
        Section { area: std::f64::consts::PI * (ro * ro - ri * ri), iy: i, iz: i, j: 2.0 * i }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Material {
    pub e: f64,
    pub g: f64,
    pub density: f64,
}

impl Material {
    pub fn isotropic(e: f64, poisson: f64, density: f64) -> Material {
        Material { e, g: e / (2.0 * (1.0 + poisson)), density }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Member {
    pub a: usize,
    pub b: usize,
    pub section: Section,
    pub material: Material,
    pub up: [f64; 3],
    pub pinned: [bool; 2],
}

pub const FIXED: [bool; 6] = [true; 6];
pub const PINNED: [bool; 6] = [true, true, true, false, false, false];

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Model {
    pub nodes: Vec<[f64; 3]>,
    pub members: Vec<Member>,
    pub supports: Vec<(usize, [bool; 6])>,
    pub masses: Vec<(usize, f64)>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Case {
    pub nodal: Vec<(usize, [f64; 6])>,
    pub uniform: Vec<(usize, [f64; 3])>,
    pub gravity: Option<[f64; 3]>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Static {
    pub displacement: Vec<[f64; 6]>,
    pub reaction: Vec<[f64; 6]>,
    pub ends: Vec<[f64; 12]>,
}

impl Static {
    pub fn combine(parts: &[(&Static, f64)]) -> Option<Static> {
        let (first, _) = parts.first()?;
        let mut out = Static {
            displacement: vec![[0.0; 6]; first.displacement.len()],
            reaction: vec![[0.0; 6]; first.reaction.len()],
            ends: vec![[0.0; 12]; first.ends.len()],
        };
        for (part, f) in parts {
            let add = |to: &mut [f64], from: &[f64]| {
                to.iter_mut().zip(from).for_each(|(t, v)| *t += f * v);
            };
            out.displacement.iter_mut().zip(&part.displacement).for_each(|(t, v)| add(t, v));
            out.reaction.iter_mut().zip(&part.reaction).for_each(|(t, v)| add(t, v));
            out.ends.iter_mut().zip(&part.ends).for_each(|(t, v)| add(t, v));
        }
        Some(out)
    }

    pub fn axial(&self, member: usize) -> f64 {
        self.ends[member][6]
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Mode {
    pub value: f64,
    pub shape: Vec<[f64; 6]>,
}

impl Mode {
    pub fn frequency(&self) -> f64 {
        self.value.sqrt() / (2.0 * std::f64::consts::PI)
    }
}

struct Element {
    local: element::Local,
    t: M12,
    k: M12,
    condensed: Condensed,
}

fn elements(model: &Model) -> Result<Vec<Element>, String> {
    let n = model.nodes.len();
    model
        .members
        .iter()
        .enumerate()
        .map(|(i, m)| {
            if m.a >= n || m.b >= n {
                return Err(format!("member {i} names a node the model does not have"));
            }
            let local = element::local(model, m).map_err(|e| format!("member {i}: {e}"))?;
            let raw = element::stiffness(m, local.length);
            let condensed = Condensed::new(&raw, &element::released(m))
                .map_err(|e| format!("member {i}: {e}"))?;
            let k = condensed.matrix(&raw);
            let t = element::transform(&local.rotation);
            Ok(Element { local, t, k, condensed })
        })
        .collect()
}

struct Dofs {
    index: Vec<usize>,
    free: usize,
}

const NONE: usize = usize::MAX;

fn number(model: &Model, diagonal: &[f64], loaded: &[bool]) -> Result<Dofs, String> {
    let total = 6 * model.nodes.len();
    let mut held = vec![false; total];
    for (node, fixed) in &model.supports {
        if *node >= model.nodes.len() {
            return Err(format!("a support names node {node}, which the model does not have"));
        }
        for (d, f) in fixed.iter().enumerate() {
            held[6 * node + d] |= *f;
        }
    }
    let scale = diagonal.iter().copied().fold(0.0f64, f64::max);
    let mut index = vec![NONE; total];
    let mut free = 0;
    for d in 0..total {
        if held[d] {
            continue;
        }
        if diagonal[d] <= 1e-12 * scale {
            if loaded[d] {
                return Err(format!("node {} is loaded in a direction nothing holds it in", d / 6));
            }
            continue;
        }
        index[d] = free;
        free += 1;
    }
    Ok(Dofs { index, free })
}

fn scatter(
    model: &Model,
    elements: &[Element],
    each: impl Fn(usize, &Element) -> M12,
) -> Vec<((usize, usize), f64)> {
    let mut out = Vec::new();
    for (i, (m, e)) in model.members.iter().zip(elements).enumerate() {
        let g = element::to_global(&e.t, &each(i, e));
        let at = |p: usize| if p < 6 { 6 * m.a + p } else { 6 * m.b + p - 6 };
        for p in 0..12 {
            for q in 0..12 {
                if g[p][q] != 0.0 {
                    out.push(((at(p), at(q)), g[p][q]));
                }
            }
        }
    }
    out
}

fn loads(model: &Model, elements: &[Element], case: &Case) -> Result<(Vec<f64>, Vec<V12>), String> {
    let total = 6 * model.nodes.len();
    let mut f = vec![0.0; total];
    let mut member: Vec<V12> = vec![[0.0; 12]; model.members.len()];
    for (node, load) in &case.nodal {
        if *node >= model.nodes.len() {
            return Err(format!("a load names node {node}, which the model does not have"));
        }
        for d in 0..6 {
            f[6 * node + d] += load[d];
        }
    }
    let mut lines: Vec<(usize, [f64; 3])> = case.uniform.clone();
    if let Some(g) = case.gravity {
        for (i, m) in model.members.iter().enumerate() {
            let w = m.material.density * m.section.area;
            lines.push((i, g.map(|v| v * w)));
        }
        for (node, mass) in &model.masses {
            for d in 0..3 {
                f[6 * node + d] += mass * g[d];
            }
        }
    }
    for (i, w) in lines {
        let e = elements
            .get(i)
            .ok_or_else(|| format!("a load names member {i}, which the model does not have"))?;
        let q = e
            .condensed
            .vector(&element::uniform(element::rotate(&e.local.rotation, w), e.local.length));
        for p in 0..12 {
            member[i][p] += q[p];
        }
    }
    for (i, (m, e)) in model.members.iter().zip(elements).enumerate() {
        let g = element::vector_to_global(&e.t, &member[i]);
        for p in 0..6 {
            f[6 * m.a + p] += g[p];
            f[6 * m.b + p] += g[p + 6];
        }
    }
    Ok((f, member))
}

fn diagonal(total: usize, entries: &[((usize, usize), f64)]) -> Vec<f64> {
    let mut d = vec![0.0; total];
    for ((i, j), v) in entries {
        if i == j {
            d[*i] += v;
        }
    }
    d
}

pub fn solve(model: &Model, case: &Case) -> Result<Static, String> {
    let _stage = sims_progress::stage("solving the frame", None, "");
    let elements = elements(model)?;
    let total = 6 * model.nodes.len();
    let k = scatter(model, &elements, |_, e| e.k);
    let (f, member) = loads(model, &elements, case)?;
    let loaded: Vec<bool> = f.iter().map(|v| *v != 0.0).collect();
    let dofs = number(model, &diagonal(total, &k), &loaded)?;
    let mut lower: std::collections::HashMap<(usize, usize), f64> = Default::default();
    for ((i, j), v) in &k {
        let (a, b) = (dofs.index[*i], dofs.index[*j]);
        if a != NONE && b != NONE && b <= a {
            *lower.entry((a, b)).or_default() += v;
        }
    }
    let triplets: Vec<faer::sparse::Triplet<usize, usize, f64>> =
        lower.into_iter().map(|((r, c), v)| faer::sparse::Triplet::new(r, c, v)).collect();
    let n = dofs.free;
    let mut x = vec![0.0; total];
    if n > 0 {
        let a = faer::sparse::SparseColMat::<usize, f64>::try_new_from_triplets(n, n, &triplets)
            .map_err(|e| format!("{e:?}"))?;
        let llt = a.sp_cholesky(faer::Side::Lower).map_err(|_| {
            "the frame is a mechanism: some part of it can move without straining a member; add supports or stiffer joints".to_string()
        })?;
        let mut b = faer::Mat::<f64>::zeros(n, 1);
        for d in 0..total {
            if dofs.index[d] != NONE {
                b[(dofs.index[d], 0)] = f[d];
            }
        }
        let u = llt.solve(&b);
        for d in 0..total {
            if dofs.index[d] != NONE {
                x[d] = u[(dofs.index[d], 0)];
            }
        }
    }
    let mut ku = vec![0.0; total];
    for ((i, j), v) in &k {
        ku[*i] += v * x[*j];
    }
    let mut reaction = vec![[0.0; 6]; model.nodes.len()];
    for (node, fixed) in &model.supports {
        for d in 0..6 {
            if fixed[d] {
                reaction[*node][d] = ku[6 * node + d] - f[6 * node + d];
            }
        }
    }
    let ends = model
        .members
        .iter()
        .zip(&elements)
        .zip(&member)
        .map(|((m, e), q)| {
            let mut ug = [0.0; 12];
            for p in 0..6 {
                ug[p] = x[6 * m.a + p];
                ug[p + 6] = x[6 * m.b + p];
            }
            let ul = element::vector_to_local(&e.t, &ug);
            std::array::from_fn(|i| (0..12).map(|j| e.k[i][j] * ul[j]).sum::<f64>() - q[i])
        })
        .collect();
    let displacement =
        (0..model.nodes.len()).map(|n| std::array::from_fn(|d| x[6 * n + d])).collect();
    Ok(Static { displacement, reaction, ends })
}

fn dense(n: usize, dofs: &Dofs, entries: &[((usize, usize), f64)]) -> faer::Mat<f64> {
    let mut m = faer::Mat::<f64>::zeros(n, n);
    for ((i, j), v) in entries {
        let (a, b) = (dofs.index[*i], dofs.index[*j]);
        if a != NONE && b != NONE {
            m[(a, b)] += v;
        }
    }
    m
}

const DENSE_LIMIT: usize = 8000;

#[cfg(not(target_arch = "wasm32"))]
fn par() -> faer::Par {
    faer::Par::rayon(0)
}

#[cfg(target_arch = "wasm32")]
fn par() -> faer::Par {
    faer::Par::Seq
}

fn smallest(
    model: &Model,
    stiffness: &[((usize, usize), f64)],
    other: &[((usize, usize), f64)],
    count: usize,
    stage: &sims_progress::Stage,
) -> Result<Vec<Mode>, String> {
    let total = 6 * model.nodes.len();
    let dofs = number(model, &diagonal(total, stiffness), &vec![false; total])?;
    let n = dofs.free;
    if n == 0 {
        return Ok(Vec::new());
    }
    if n > DENSE_LIMIT {
        return Err(format!(
            "{n} free degrees of freedom is more than the dense eigensolver takes ({DENSE_LIMIT})"
        ));
    }
    let k = dense(n, &dofs, stiffness);
    let b = dense(n, &dofs, other);
    stage.detail(|| format!("{n} degrees of freedom"));
    let llt = k.llt(faer::Side::Lower).map_err(|_| {
        "the frame is a mechanism: some part of it can move without straining a member".to_string()
    })?;
    let l = llt.L();
    let mut c = b;
    faer::linalg::triangular_solve::solve_lower_triangular_in_place(l, c.as_mut(), par());
    let mut c = c.transpose().to_owned();
    faer::linalg::triangular_solve::solve_lower_triangular_in_place(l, c.as_mut(), par());
    let sym = faer::Mat::<f64>::from_fn(n, n, |i, j| 0.5 * (c[(i, j)] + c[(j, i)]));
    let eig = sym.self_adjoint_eigen(faer::Side::Lower).map_err(|e| format!("{e:?}"))?;
    let (s, u) = (eig.S(), eig.U());
    let mut out = Vec::new();
    for col in (0..n).rev() {
        let mu = s[col];
        if mu <= 0.0 || out.len() >= count {
            break;
        }
        let mut y = faer::Mat::<f64>::from_fn(n, 1, |i, _| u[(i, col)]);
        faer::linalg::triangular_solve::solve_upper_triangular_in_place(
            l.transpose(),
            y.as_mut(),
            faer::Par::Seq,
        );
        let peak = (0..n).map(|i| y[(i, 0)].abs()).fold(0.0f64, f64::max).max(1e-300);
        let mut shape = vec![[0.0; 6]; model.nodes.len()];
        for d in 0..total {
            if dofs.index[d] != NONE {
                shape[d / 6][d % 6] = y[(dofs.index[d], 0)] / peak;
            }
        }
        out.push(Mode { value: 1.0 / mu, shape });
    }
    Ok(out)
}

pub fn modes(model: &Model, count: usize) -> Result<Vec<Mode>, String> {
    let stage = sims_progress::stage("finding the frame's modes", None, "");
    let elements = elements(model)?;
    let k = scatter(model, &elements, |_, e| e.k);
    let mut m = scatter(model, &elements, |i, e| {
        e.condensed.matrix(&element::mass(&model.members[i], e.local.length))
    });
    for (node, mass) in &model.masses {
        for d in 0..3 {
            m.push(((6 * node + d, 6 * node + d), *mass));
        }
    }
    smallest(model, &k, &m, count, &stage)
}

pub fn buckling(model: &Model, case: &Case, count: usize) -> Result<Vec<Mode>, String> {
    let stage = sims_progress::stage("finding the frame's buckling loads", None, "");
    let pre = solve(model, case)?;
    let elements = elements(model)?;
    let k = scatter(model, &elements, |_, e| e.k);
    let kg = scatter(model, &elements, |i, e| {
        let g = element::geometric(-pre.axial(i), e.local.length);
        e.condensed.matrix(&g)
    });
    smallest(model, &k, &kg, count, &stage)
}
