#![cfg(not(target_arch = "wasm32"))]

use fieldgrid::gpu::{Gpu, Options, Software};

const DOUBLE: &str = "
@group(0) @binding(0) var<storage, read_write> data: array<f32>;

@compute @workgroup_size(64)
fn double(@builtin(global_invocation_id) id: vec3<u32>) {
    if id.x < arrayLength(&data) {
        data[id.x] = data[id.x] * 2.0;
    }
}
";

fn open() -> Option<Gpu> {
    let options = Options { software: Software::Allow, ..Options::default() };
    match Gpu::open_blocking(&options) {
        Ok(g) => Some(g),
        Err(e) if !cfg!(feature = "gpu-tests") => {
            eprintln!("skipped: {e}");
            None
        }
        Err(e) => panic!("{e}"),
    }
}

#[test]
fn a_kernel_runs_and_reads_back() {
    let Some(gpu) = open() else {
        return;
    };
    let input: Vec<f32> = (0..1000).map(|i| i as f32).collect();
    let buffer = gpu.storage("data", &input);
    let pipeline = gpu.pipeline(DOUBLE, "double");
    let bind = gpu.bind(&pipeline, &[&buffer]);
    let mut encoder = gpu.device.create_command_encoder(&Default::default());
    {
        let mut pass = encoder.begin_compute_pass(&Default::default());
        pass.set_pipeline(&pipeline);
        pass.set_bind_group(0, &bind, &[]);
        pass.dispatch_workgroups(1000u32.div_ceil(64), 1, 1);
    }
    let out: Vec<f32> = gpu.take_now(gpu.fetch(encoder, &buffer, 0..4000)).unwrap();
    assert_eq!(out, input.iter().map(|x| x * 2.0).collect::<Vec<_>>());

    let tail: Vec<f32> = gpu.read(&buffer, 990..1000).unwrap();
    assert_eq!(tail, (990..1000).map(|i| 2.0 * i as f32).collect::<Vec<_>>());

    let encoder = gpu.device.create_command_encoder(&Default::default());
    let before: Vec<f32> = gpu.take_now(gpu.fetch_cleared(encoder, &buffer)).unwrap();
    assert_eq!(before[999], 1998.0);
    let after: Vec<f32> = gpu.read_all(&buffer).unwrap();
    assert!(after.iter().all(|x| *x == 0.0));
}

#[test]
fn a_uniform_and_a_zeroed_buffer_have_their_sizes() {
    let Some(gpu) = open() else {
        return;
    };
    assert_eq!(gpu.zeroed("z", 0).size(), 4);
    assert_eq!(gpu.zeroed("z", 400).size(), 400);
    assert_eq!(gpu.uniform("u", &[1u32, 2, 3, 4]).size(), 16);
    assert!(gpu.max_binding() > 0);
    assert!(!gpu.name.is_empty());
}
