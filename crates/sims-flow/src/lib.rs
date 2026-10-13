pub mod bench;
pub mod lattice;
pub mod rotor;
pub mod voxel;

use anyhow::{Result, bail};
use sims_field::Gpu;
use sims_field::gpu::{Options, Reading, Software};

pub mod geom {
    pub type Point3 = cgmath::Point3<f64>;
    pub type Vector3 = cgmath::Vector3<f64>;
    pub use cgmath::prelude::*;
}

pub fn open(prefer_cpu: bool) -> Result<Gpu> {
    if cfg!(target_arch = "wasm32") {
        bail!("flow simulations do not run in the browser");
    }
    let options = Options {
        label: "flow",
        backends: wgpu::Backends::all(),
        software: if prefer_cpu { Software::Prefer } else { Software::Allow },
        storage_buffers: 40,
        features: wgpu::Features::SHADER_F16,
    };
    open_blocking(&options)
}

#[cfg(not(target_arch = "wasm32"))]
fn open_blocking(options: &Options) -> Result<Gpu> {
    Gpu::open_blocking(options).map_err(anyhow::Error::msg)
}

#[cfg(target_arch = "wasm32")]
fn open_blocking(_: &Options) -> Result<Gpu> {
    bail!("flow simulations do not run in the browser")
}

fn fetch(gpu: &Gpu, encoder: wgpu::CommandEncoder, buffer: &wgpu::Buffer) -> Reading {
    gpu.fetch(encoder, buffer, 0..buffer.size())
}

#[cfg(not(target_arch = "wasm32"))]
fn read(gpu: &Gpu, buffer: &wgpu::Buffer) -> Result<Vec<f32>> {
    gpu.read_all(buffer).map_err(anyhow::Error::msg)
}

#[cfg(not(target_arch = "wasm32"))]
fn take(gpu: &Gpu, reading: Reading) -> Result<Vec<f32>> {
    gpu.take_now(reading).map_err(anyhow::Error::msg)
}

#[cfg(target_arch = "wasm32")]
fn read(_: &Gpu, _: &wgpu::Buffer) -> Result<Vec<f32>> {
    bail!("flow simulations do not run in the browser")
}

#[cfg(target_arch = "wasm32")]
fn take(_: &Gpu, _: Reading) -> Result<Vec<f32>> {
    bail!("flow simulations do not run in the browser")
}

#[cfg(not(target_arch = "wasm32"))]
fn wait(gpu: &Gpu, index: Option<wgpu::SubmissionIndex>) -> Result<()> {
    gpu.wait(index).map_err(anyhow::Error::msg)
}

#[cfg(target_arch = "wasm32")]
fn wait(_: &Gpu, _: Option<wgpu::SubmissionIndex>) -> Result<()> {
    bail!("flow simulations do not run in the browser")
}
