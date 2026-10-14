use crate::aggregate::Rigid;
use crate::solve::{Block, Blocks, Solution, TwoLevel};
use bytemuck::{Pod, Zeroable};
use sims_field::Gpu;
use sims_field::gpu::{Options, Software};

pub fn open(software: Software) -> Result<Gpu, String> {
    let options = Options {
        label: "sims-fem",
        software,
        features: wgpu::Features::SHADER_F64,
        ..Default::default()
    };
    let gpu = Gpu::open_blocking(&options)?;
    if !gpu.device.features().contains(wgpu::Features::SHADER_F64) {
        return Err(format!("{} has no 64-bit floats in shaders", gpu.name));
    }
    Ok(gpu)
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Range {
    first: u32,
    count: u32,
    nodes: u32,
    pad: u32,
}

const FINE_WEIGHT: usize = 0;
const ALPHA: usize = 1;
const LESS_ALPHA: usize = 2;
const BETA: usize = 3;
const COARSE_WEIGHT: usize = 4;
const PARTIALS: u32 = 1024;

const KERNELS: [&str; 12] = [
    "multiply",
    "project",
    "sweep",
    "jacobi",
    "clusters",
    "gather",
    "prolong",
    "add_scaled",
    "scale_add",
    "dot",
    "aggregate",
    "spread",
];

fn groups(threads: usize) -> u32 {
    (threads as u32).div_ceil(256).max(1)
}

fn floats(blocks: &[Block]) -> Vec<f64> {
    blocks.iter().flatten().copied().collect()
}

fn nonempty<T: Clone>(v: Vec<T>, fill: T) -> Vec<T> {
    if v.is_empty() { vec![fill; 4] } else { v }
}

struct Step {
    pipeline: usize,
    group: wgpu::BindGroup,
    size: u32,
}

enum Op<'s> {
    Run(&'s Step),
    Clear(&'s wgpu::Buffer),
    Copy(&'s wgpu::Buffer, &'s wgpu::Buffer),
}

struct Part {
    first: usize,
    count: usize,
    starts: wgpu::Buffer,
    columns: wgpu::Buffer,
    values: wgpu::Buffer,
}

const SLICE: usize = 32;

fn partitions(gpu: &Gpu, a: &Blocks, bytes: u64) -> Vec<Part> {
    let limit = (bytes / (9 * 8)).max(1) as usize;
    let nodes = a.nodes();
    let width = |s: usize| {
        (s * SLICE..((s + 1) * SLICE).min(nodes))
            .map(|r| a.starts[r + 1] - a.starts[r])
            .max()
            .unwrap_or(0)
    };
    let slices = nodes.div_ceil(SLICE);
    let mut out = Vec::new();
    let mut first = 0;
    while first < slices {
        let mut end = first;
        let mut size = 0;
        while end < slices && (end == first || size + width(end) * SLICE <= limit) {
            size += width(end) * SLICE;
            end += 1;
        }
        let mut info = Vec::with_capacity(2 * (end - first));
        let mut columns = Vec::with_capacity(size);
        let mut values: Vec<f64> = Vec::with_capacity(9 * size);
        for slice in first..end {
            let w = width(slice);
            info.extend_from_slice(&[columns.len() as u32, w as u32]);
            for k in 0..w {
                let entries: Vec<Option<usize>> = (0..SLICE)
                    .map(|lane| {
                        let row = slice * SLICE + lane;
                        (row < nodes).then(|| a.starts[row] + k).filter(|&e| e < a.starts[row + 1])
                    })
                    .collect();
                for (lane, at) in entries.iter().enumerate() {
                    let row = (slice * SLICE + lane).min(nodes - 1);
                    columns.push(at.map_or(row as u32, |e| a.columns[e]));
                }
                for j in 0..9 {
                    values.extend(entries.iter().map(|at| at.map_or(0.0, |e| a.values[e][j])));
                }
            }
        }
        let row_start = first * SLICE;
        let row_end = (end * SLICE).min(nodes);
        out.push(Part {
            first: row_start,
            count: row_end - row_start,
            starts: gpu.storage("slices", &info),
            columns: gpu.storage("columns", &nonempty(columns, 0)),
            values: gpu.storage("values", &nonempty(values, 0.0)),
        });
        first = end;
    }
    out
}

struct Builder<'g> {
    gpu: &'g Gpu,
    pipelines: &'g [wgpu::ComputePipeline],
}

impl Builder<'_> {
    fn range(&self, first: usize, count: usize, nodes: usize) -> wgpu::Buffer {
        let r = Range { first: first as u32, count: count as u32, nodes: nodes as u32, pad: 0 };
        self.gpu.uniform("range", &r)
    }

    fn step(&self, name: &str, threads: usize, entries: &[(u32, &wgpu::Buffer)]) -> Step {
        let pipeline = KERNELS.iter().position(|k| *k == name).expect("a kernel");
        let entries: Vec<wgpu::BindGroupEntry> = entries
            .iter()
            .map(|(binding, b)| wgpu::BindGroupEntry {
                binding: *binding,
                resource: b.as_entire_binding(),
            })
            .collect();
        let group = self.gpu.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some(name),
            layout: &self.pipelines[pipeline].get_bind_group_layout(0),
            entries: &entries,
        });
        Step { pipeline, group, size: groups(threads) }
    }

    fn multiply(
        &self,
        parts: &[Part],
        nodes: usize,
        x: &wgpu::Buffer,
        y: &wgpu::Buffer,
    ) -> Vec<Step> {
        parts
            .iter()
            .map(|part| {
                let r = self.range(part.first, part.count, nodes);
                self.step(
                    "multiply",
                    part.count,
                    &[
                        (0, &r),
                        (1, &part.starts),
                        (2, &part.columns),
                        (3, &part.values),
                        (4, x),
                        (5, y),
                    ],
                )
            })
            .collect()
    }
}

struct Vectors {
    x: wgpu::Buffer,
    b: wgpu::Buffer,
    ax: wgpu::Buffer,
    res: wgpu::Buffer,
}

impl Vectors {
    fn new(gpu: &Gpu, nodes: usize) -> Vectors {
        let bytes = (3 * nodes.max(1)) as u64 * 8;
        Vectors {
            x: gpu.zeroed("x", bytes),
            b: gpu.zeroed("b", bytes),
            ax: gpu.zeroed("ax", bytes),
            res: gpu.zeroed("res", bytes),
        }
    }
}

struct Level {
    smooth: Vec<Step>,
    residual: Vec<Step>,
    project_x: Step,
    project_b: Step,
}

struct Matrix<'m> {
    parts: &'m [Part],
    nodes: usize,
    projectors: &'m wgpu::Buffer,
    inverse: &'m wgpu::Buffer,
}

fn level(
    build: &Builder,
    m: &Matrix,
    weight: usize,
    v: &Vectors,
    scalars: &wgpu::Buffer,
    clusters: Option<Step>,
) -> Level {
    let whole = build.range(weight, m.nodes, m.nodes);
    let sweep = || {
        build.step(
            "sweep",
            m.nodes,
            &[(0, &whole), (3, m.projectors), (4, &v.ax), (6, &v.b), (5, &v.res)],
        )
    };
    let project = |t: &wgpu::Buffer| {
        build.step("project", m.nodes, &[(0, &whole), (3, m.projectors), (5, t)])
    };
    let mut residual = build.multiply(m.parts, m.nodes, &v.x, &v.ax);
    residual.push(sweep());
    let mut smooth = build.multiply(m.parts, m.nodes, &v.x, &v.ax);
    smooth.push(sweep());
    smooth.push(build.step(
        "jacobi",
        m.nodes,
        &[(0, &whole), (3, m.inverse), (4, &v.res), (5, &v.x), (7, scalars)],
    ));
    smooth.extend(clusters);
    smooth.push(project(&v.x));
    Level { smooth, residual, project_x: project(&v.x), project_b: project(&v.b) }
}

struct Aggregates {
    gather: Step,
    spread: Step,
    sums: wgpu::Buffer,
    solved: wgpu::Buffer,
    count: usize,
}

fn aggregates_of(build: &Builder, rigid: &Rigid, v: &Vectors) -> Aggregates {
    let gpu = build.gpu;
    let nodes = rigid.a.nodes();
    let mut members: Vec<Vec<u32>> = vec![Vec::new(); rigid.count];
    for (node, &g) in rigid.owner.iter().enumerate() {
        members[g as usize].push(node as u32);
    }
    let mut starts = vec![0u32];
    let mut list = Vec::new();
    for m in members {
        list.extend(m);
        starts.push(list.len() as u32);
    }
    let starts = gpu.storage("aggregate starts", &starts);
    let list = gpu.storage("aggregate members", &nonempty(list, 0));
    let owner = gpu.storage("owner", &nonempty(rigid.owner.clone(), 0));
    let modes: Vec<f64> = rigid.modes.iter().flatten().copied().collect();
    let modes = gpu.storage("modes", &nonempty(modes, 0.0));
    let sums = gpu.zeroed("aggregate sums", (6 * rigid.count.max(1)) as u64 * 8);
    let solved = gpu.zeroed("aggregate solution", (6 * rigid.count.max(1)) as u64 * 8);
    let range = build.range(0, rigid.count, nodes);
    let gather = build.step(
        "aggregate",
        6 * rigid.count,
        &[(0, &range), (1, &starts), (2, &list), (3, &modes), (4, &v.res), (5, &sums)],
    );
    let spread = build.step(
        "spread",
        nodes,
        &[(0, &range), (1, &owner), (3, &modes), (4, &solved), (5, &v.x)],
    );
    Aggregates { gather, spread, sums, solved, count: rigid.count }
}

fn children_of(level: &TwoLevel) -> (Vec<u32>, Vec<u32>, Vec<f64>) {
    let mut lists: Vec<Vec<(u32, f64)>> = vec![Vec::new(); level.coarse.nodes];
    for (i, parents) in level.coarse.parents.iter().enumerate() {
        for &(p, w) in parents {
            if w != 0.0 {
                lists[p as usize].push((i as u32, w));
            }
        }
    }
    let mut starts = vec![0u32];
    let mut children = Vec::new();
    let mut weights = Vec::new();
    for list in lists {
        for (c, w) in list {
            children.push(c);
            weights.push(w);
        }
        starts.push(children.len() as u32);
    }
    (starts, children, weights)
}

fn cluster_step(
    build: &Builder,
    level: &TwoLevel,
    v: &Vectors,
    scalars: &wgpu::Buffer,
) -> Option<Step> {
    if level.clusters.is_empty() {
        return None;
    }
    let mut meta = Vec::new();
    let mut members = Vec::new();
    let mut mats = Vec::new();
    for cluster in &level.clusters {
        let size = 3 * cluster.nodes.len();
        let first_node = members.len() as u32;
        members.extend_from_slice(&cluster.nodes);
        for row in 0..size {
            meta.extend_from_slice(&[
                first_node,
                (mats.len() + row * size) as u32,
                size as u32,
                row as u32,
            ]);
        }
        mats.extend_from_slice(&cluster.inverse);
    }
    let rows = meta.len() / 4;
    let gpu = build.gpu;
    let meta = gpu.storage("cluster rows", &meta);
    let members = gpu.storage("cluster nodes", &members);
    let mats = gpu.storage("cluster inverses", &mats);
    let range = build.range(FINE_WEIGHT, rows, level.a.nodes());
    Some(build.step(
        "clusters",
        rows,
        &[(0, &range), (1, &meta), (2, &members), (3, &mats), (4, &v.res), (5, &v.x), (7, scalars)],
    ))
}

pub struct Solver<'a> {
    gpu: &'a Gpu,
    level: TwoLevel<'a>,
    pipelines: Vec<wgpu::ComputePipeline>,
    nodes: usize,
    scalars: wgpu::Buffer,
    fine: Vectors,
    coarse: Vectors,
    x: wgpu::Buffer,
    r: wgpu::Buffer,
    p: wgpu::Buffer,
    partials: wgpu::Buffer,
    fine_level: Level,
    coarse_level: Option<Level>,
    aggregates: Option<Aggregates>,
    gather: Step,
    prolong: Step,
    multiply_p: Vec<Step>,
    project_q: Step,
    add_x: Step,
    add_r: Step,
    update_p: Step,
    dot_rz: Step,
    dot_pq: Step,
    dot_rr: Step,
}

impl<'a> Solver<'a> {
    pub fn new(gpu: &'a Gpu, level: TwoLevel<'a>) -> Result<Solver<'a>, String> {
        Solver::with_partition_bytes(gpu, level, gpu.max_binding())
    }

    pub fn with_partition_bytes(
        gpu: &'a Gpu,
        level: TwoLevel<'a>,
        bytes: u64,
    ) -> Result<Solver<'a>, String> {
        let _stage = sims_progress::stage("moving the stiffness to the GPU", None, "");
        let module = gpu.device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("sims-fem"),
            source: wgpu::ShaderSource::Wgsl(include_str!("gpu.wgsl").into()),
        });
        let pipelines: Vec<wgpu::ComputePipeline> = KERNELS
            .iter()
            .map(|entry| {
                gpu.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some(entry),
                    layout: None,
                    module: &module,
                    entry_point: Some(entry),
                    compilation_options: Default::default(),
                    cache: None,
                })
            })
            .collect();
        let build = Builder { gpu, pipelines: &pipelines };
        let nodes = level.a.nodes();
        let coarse_nodes = level.coarse.nodes;
        let rigid_weight = level.rigid.as_ref().map_or(0.0, |r| r.weight);
        let scalars =
            gpu.storage("scalars", &[level.weight, 0.0, 0.0, 0.0, rigid_weight, 0.0, 0.0, 0.0]);
        let fine = Vectors::new(gpu, nodes);
        let coarse = Vectors::new(gpu, coarse_nodes);
        let bytes3 = (3 * nodes.max(1)) as u64 * 8;
        let (x, r, p, q) = (
            gpu.zeroed("x", bytes3),
            gpu.zeroed("r", bytes3),
            gpu.zeroed("p", bytes3),
            gpu.zeroed("q", bytes3),
        );
        let partials = gpu.zeroed("partials", PARTIALS as u64 * 8);
        let parts = partitions(gpu, level.a, bytes);
        let projectors = gpu.storage("projectors", &nonempty(floats(level.projectors), 0.0));
        let inverse = gpu.storage("inverse", &nonempty(floats(&level.inverse), 0.0));
        let fine_matrix =
            Matrix { parts: &parts, nodes, projectors: &projectors, inverse: &inverse };
        let clusters = cluster_step(&build, &level, &fine, &scalars);
        let fine_level = self::level(&build, &fine_matrix, FINE_WEIGHT, &fine, &scalars, clusters);
        let (coarse_level, aggregates) = match &level.rigid {
            Some(rigid) => {
                let rigid_parts = partitions(gpu, &rigid.a, bytes);
                let rp = gpu.storage("projectors", &nonempty(floats(&rigid.projectors), 0.0));
                let ri = gpu.storage("inverse", &nonempty(floats(&rigid.inverse), 0.0));
                let m = Matrix {
                    parts: &rigid_parts,
                    nodes: rigid.a.nodes(),
                    projectors: &rp,
                    inverse: &ri,
                };
                let l = self::level(&build, &m, COARSE_WEIGHT, &coarse, &scalars, None);
                (Some(l), Some(aggregates_of(&build, rigid, &coarse)))
            }
            None => (None, None),
        };
        let (child_starts, children, child_weights) = children_of(&level);
        let child_starts = gpu.storage("child starts", &child_starts);
        let children = gpu.storage("children", &nonempty(children, 0));
        let child_weights = gpu.storage("child weights", &nonempty(child_weights, 0.0));
        let coarse_range = build.range(0, coarse_nodes, coarse_nodes);
        let gather = build.step(
            "gather",
            coarse_nodes,
            &[
                (0, &coarse_range),
                (1, &child_starts),
                (2, &children),
                (3, &child_weights),
                (4, &fine.res),
                (5, &coarse.b),
            ],
        );
        let parent_nodes: Vec<u32> =
            level.coarse.parents.iter().flat_map(|p| [p[0].0, p[1].0]).collect();
        let parent_weights: Vec<f64> =
            level.coarse.parents.iter().flat_map(|p| [p[0].1, p[1].1]).collect();
        let parent_nodes = gpu.storage("parents", &nonempty(parent_nodes, 0));
        let parent_weights = gpu.storage("parent weights", &nonempty(parent_weights, 0.0));
        let whole = build.range(0, nodes, nodes);
        let prolong = build.step(
            "prolong",
            nodes,
            &[(0, &whole), (1, &parent_nodes), (3, &parent_weights), (4, &coarse.x), (5, &fine.x)],
        );
        let multiply_p = build.multiply(&parts, nodes, &p, &q);
        let project_q = build.step("project", nodes, &[(0, &whole), (3, &projectors), (5, &q)]);
        let scaled = |first: usize, src: &wgpu::Buffer, dst: &wgpu::Buffer, name: &str| {
            let rr = build.range(first, 0, nodes);
            build.step(name, 3 * nodes, &[(0, &rr), (4, src), (5, dst), (7, &scalars)])
        };
        let dots = build.range(0, PARTIALS as usize, nodes);
        let dot = |a: &wgpu::Buffer, b: &wgpu::Buffer| {
            let mut step = build.step("dot", 0, &[(0, &dots), (4, a), (6, b), (5, &partials)]);
            step.size = PARTIALS;
            step
        };
        let add_x = scaled(ALPHA, &p, &x, "add_scaled");
        let add_r = scaled(LESS_ALPHA, &q, &r, "add_scaled");
        let update_p = scaled(BETA, &fine.x, &p, "scale_add");
        let (dot_rz, dot_pq, dot_rr) = (dot(&r, &fine.x), dot(&p, &q), dot(&r, &r));
        Ok(Solver {
            gpu,
            level,
            pipelines,
            nodes,
            scalars,
            fine,
            coarse,
            x,
            r,
            p,
            partials,
            fine_level,
            coarse_level,
            aggregates,
            gather,
            prolong,
            multiply_p,
            project_q,
            add_x,
            add_r,
            update_p,
            dot_rz,
            dot_pq,
            dot_rr,
        })
    }

    fn encode(&self, ops: &[Op]) -> wgpu::CommandEncoder {
        let mut encoder = self.gpu.device.create_command_encoder(&Default::default());
        let mut i = 0;
        while i < ops.len() {
            match ops[i] {
                Op::Clear(b) => {
                    encoder.clear_buffer(b, 0, None);
                    i += 1;
                }
                Op::Copy(from, to) => {
                    encoder.copy_buffer_to_buffer(from, 0, to, 0, from.size().min(to.size()));
                    i += 1;
                }
                Op::Run(_) => {
                    let mut pass = encoder.begin_compute_pass(&Default::default());
                    while let Some(Op::Run(step)) = ops.get(i) {
                        pass.set_pipeline(&self.pipelines[step.pipeline]);
                        pass.set_bind_group(0, &step.group, &[]);
                        pass.dispatch_workgroups(
                            step.size.min(65535),
                            step.size.div_ceil(65535),
                            1,
                        );
                        i += 1;
                    }
                }
            }
        }
        encoder
    }

    fn submit(&self, ops: &[Op]) {
        self.gpu.queue.submit([self.encode(ops).finish()]);
    }

    fn read(&self, ops: &[Op], buffer: &wgpu::Buffer, count: usize) -> Result<Vec<f64>, String> {
        let reading = self.gpu.fetch(self.encode(ops), buffer, 0..count as u64 * 8);
        self.gpu.take_now(reading)
    }

    fn dot<'s>(&'s self, mut ops: Vec<Op<'s>>, step: &'s Step) -> Result<f64, String> {
        ops.push(Op::Run(step));
        Ok(self.read(&ops, &self.partials, PARTIALS as usize)?.iter().sum())
    }

    fn set(&self, index: usize, value: f64) {
        self.gpu.queue.write_buffer(&self.scalars, index as u64 * 8, bytemuck::bytes_of(&value));
    }

    fn smooth_twice<'s>(level: &'s Level, ops: &mut Vec<Op<'s>>) {
        for _ in 0..2 {
            ops.extend(level.smooth.iter().map(Op::Run));
        }
    }

    fn precondition(&self) -> Result<(), String> {
        let (f, c) = (&self.fine_level, &self.coarse_level);
        let mut ops =
            vec![Op::Copy(&self.r, &self.fine.b), Op::Clear(&self.fine.x), Op::Run(&f.project_b)];
        Self::smooth_twice(f, &mut ops);
        ops.extend(f.residual.iter().map(Op::Run));
        ops.push(Op::Run(&self.gather));
        match (c, &self.aggregates, self.level.rigid.as_deref()) {
            (Some(c), Some(agg), Some(rigid)) => {
                ops.push(Op::Run(&c.project_b));
                ops.push(Op::Clear(&self.coarse.x));
                Self::smooth_twice(c, &mut ops);
                ops.extend(c.residual.iter().map(Op::Run));
                ops.push(Op::Run(&agg.gather));
                let sums = self.read(&ops, &agg.sums, 6 * agg.count)?;
                let solved = rigid.solve_aggregates(&sums);
                if !solved.is_empty() {
                    self.gpu.queue.write_buffer(&agg.solved, 0, bytemuck::cast_slice(&solved));
                }
                let mut ops = vec![Op::Run(&agg.spread), Op::Run(&c.project_x)];
                Self::smooth_twice(c, &mut ops);
                self.finish(ops);
            }
            _ => {
                let gathered = self.read(&ops, &self.coarse.b, 3 * self.level.coarse.nodes)?;
                let correction = self.level.coarse_correction(&gathered);
                if !correction.is_empty() {
                    self.gpu.queue.write_buffer(
                        &self.coarse.x,
                        0,
                        bytemuck::cast_slice(&correction),
                    );
                }
                self.finish(Vec::new());
            }
        }
        Ok(())
    }

    fn finish<'s>(&'s self, mut ops: Vec<Op<'s>>) {
        ops.push(Op::Run(&self.prolong));
        ops.push(Op::Run(&self.fine_level.project_x));
        Self::smooth_twice(&self.fine_level, &mut ops);
        self.submit(&ops);
    }

    pub fn precondition_vector(&self, r: &[f64]) -> Result<Vec<f64>, String> {
        self.gpu.queue.write_buffer(&self.r, 0, bytemuck::cast_slice(r));
        self.precondition()?;
        self.gpu.read(&self.fine.x, 0..3 * self.nodes)
    }

    pub fn level(&self) -> &TwoLevel<'a> {
        &self.level
    }

    pub fn aggregated(&self) -> bool {
        self.aggregates.is_some()
    }

    pub fn solve(&self, b: &[f64], tolerance: f64, limit: usize) -> Result<Solution, String> {
        let n = 3 * self.nodes;
        let mut rhs = b.to_vec();
        crate::solve::project(self.level.projectors, &mut rhs);
        let norm = crate::solve::dot(&rhs, &rhs).sqrt();
        if norm == 0.0 {
            return Ok(Solution { x: vec![0.0; n], iterations: 0, residual: 0.0 });
        }
        self.gpu.queue.write_buffer(&self.r, 0, bytemuck::cast_slice(&rhs));
        self.submit(&[Op::Clear(&self.x)]);
        self.precondition()?;
        let mut rz = self.dot(vec![Op::Copy(&self.fine.x, &self.p)], &self.dot_rz)?;
        let mut iterations = 0;
        let mut residual = 1.0;
        let mut best = 1.0f64;
        let stage = sims_progress::stage("solving on the GPU", None, "iterations");
        while iterations < limit {
            let mut ops: Vec<Op> = self.multiply_p.iter().map(Op::Run).collect();
            ops.push(Op::Run(&self.project_q));
            let pq = self.dot(ops, &self.dot_pq)?;
            if pq <= 0.0 {
                break;
            }
            let alpha = rz / pq;
            self.set(ALPHA, alpha);
            self.set(LESS_ALPHA, -alpha);
            let rr = self.dot(vec![Op::Run(&self.add_x), Op::Run(&self.add_r)], &self.dot_rr)?;
            iterations += 1;
            residual = rr.max(0.0).sqrt() / norm;
            best = best.min(residual);
            stage.fraction(sims_progress::residual_fraction(1.0, best, tolerance));
            stage.progress(iterations as f64, || {
                format!("residual {residual:.1e}, aiming for {tolerance:.0e}")
            });
            if residual < tolerance {
                break;
            }
            self.precondition()?;
            let next = self.dot(Vec::new(), &self.dot_rz)?;
            self.set(BETA, next / rz);
            rz = next;
            self.submit(&[Op::Run(&self.update_p)]);
        }
        let x: Vec<f64> = self.gpu.read(&self.x, 0..n)?;
        Ok(Solution { x, iterations, residual })
    }
}
