use std::future::Future;
use std::ops::Range;
use std::pin::Pin;
use std::sync::{Arc, Mutex};
use std::task::{Context, Poll, Waker};
use wgpu::util::DeviceExt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Software {
    Refuse,
    Allow,
    Prefer,
}

#[derive(Clone, Debug)]
pub struct Options {
    pub label: &'static str,
    pub backends: wgpu::Backends,
    pub software: Software,
    pub storage_buffers: u32,
    pub features: wgpu::Features,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            label: "sims-field",
            backends: wgpu::Backends::PRIMARY,
            software: Software::Refuse,
            storage_buffers: 8,
            features: wgpu::Features::empty(),
        }
    }
}

pub struct Gpu {
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
    pub info: wgpu::AdapterInfo,
    pub name: String,
}

pub struct Reading {
    staging: wgpu::Buffer,
    index: wgpu::SubmissionIndex,
}

fn rank(software: Software, info: &wgpu::AdapterInfo) -> (usize, usize) {
    let kind = match info.device_type {
        wgpu::DeviceType::DiscreteGpu => 0,
        wgpu::DeviceType::IntegratedGpu => 1,
        wgpu::DeviceType::Cpu => 2,
        _ => 3,
    };
    let kind = if software == Software::Prefer { [2, 1, 0, 3][kind] } else { kind };
    let backend = usize::from(info.backend != wgpu::Backend::Vulkan);
    (kind, backend)
}

pub async fn adapter(options: &Options) -> Result<wgpu::Adapter, String> {
    let backends = wgpu::Backends::from_env().unwrap_or(options.backends);
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
        backends,
        ..wgpu::InstanceDescriptor::new_without_display_handle_from_env()
    });
    instance
        .enumerate_adapters(backends)
        .await
        .into_iter()
        .filter(|a| {
            options.software != Software::Refuse
                || a.get_info().device_type != wgpu::DeviceType::Cpu
        })
        .filter(|a| a.limits().max_storage_buffers_per_shader_stage >= options.storage_buffers)
        .min_by_key(|a| rank(options.software, &a.get_info()))
        .ok_or_else(|| {
            format!(
                "no {} with {} storage buffers per shader",
                if options.software == Software::Refuse { "GPU" } else { "compute device" },
                options.storage_buffers
            )
        })
}

impl Gpu {
    pub async fn open(options: &Options) -> Result<Gpu, String> {
        Gpu::on(adapter(options).await?, options).await
    }

    pub async fn on(adapter: wgpu::Adapter, options: &Options) -> Result<Gpu, String> {
        let info = adapter.get_info();
        let name = format!("{} ({:?})", info.name, info.backend);
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor {
                label: Some(options.label),
                required_features: adapter.features() & options.features,
                required_limits: adapter.limits(),
                ..Default::default()
            })
            .await
            .map_err(|e| format!("opening {name}: {e}"))?;
        Ok(Gpu { device, queue, info, name })
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn open_blocking(options: &Options) -> Result<Gpu, String> {
        pollster::block_on(Gpu::open(options))
    }

    pub fn max_binding(&self) -> u64 {
        self.device.limits().max_storage_buffer_binding_size
    }

    pub fn storage<T: bytemuck::Pod>(&self, label: &str, data: &[T]) -> wgpu::Buffer {
        self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some(label),
            contents: bytemuck::cast_slice(data),
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC
                | wgpu::BufferUsages::COPY_DST,
        })
    }

    pub fn zeroed(&self, label: &str, bytes: u64) -> wgpu::Buffer {
        self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size: bytes.max(4),
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC
                | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        })
    }

    pub fn uniform<T: bytemuck::Pod>(&self, label: &str, value: &T) -> wgpu::Buffer {
        self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some(label),
            contents: bytemuck::bytes_of(value),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        })
    }

    pub fn pipeline(&self, source: &str, entry: &str) -> wgpu::ComputePipeline {
        let module = self.device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some(entry),
            source: wgpu::ShaderSource::Wgsl(source.into()),
        });
        self.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(entry),
            layout: None,
            module: &module,
            entry_point: Some(entry),
            compilation_options: Default::default(),
            cache: None,
        })
    }

    pub fn bind(
        &self,
        pipeline: &wgpu::ComputePipeline,
        buffers: &[&wgpu::Buffer],
    ) -> wgpu::BindGroup {
        let entries: Vec<wgpu::BindGroupEntry> = buffers
            .iter()
            .enumerate()
            .map(|(i, b)| wgpu::BindGroupEntry {
                binding: i as u32,
                resource: b.as_entire_binding(),
            })
            .collect();
        self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &pipeline.get_bind_group_layout(0),
            entries: &entries,
        })
    }

    pub fn fetch(
        &self,
        mut encoder: wgpu::CommandEncoder,
        buffer: &wgpu::Buffer,
        bytes: Range<u64>,
    ) -> Reading {
        let staging = self.staging(&mut encoder, buffer, bytes);
        let index = self.queue.submit([encoder.finish()]);
        Reading { staging, index }
    }

    pub fn fetch_cleared(
        &self,
        mut encoder: wgpu::CommandEncoder,
        buffer: &wgpu::Buffer,
    ) -> Reading {
        let staging = self.staging(&mut encoder, buffer, 0..buffer.size());
        encoder.clear_buffer(buffer, 0, None);
        let index = self.queue.submit([encoder.finish()]);
        Reading { staging, index }
    }

    fn staging(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        buffer: &wgpu::Buffer,
        bytes: Range<u64>,
    ) -> wgpu::Buffer {
        let size = bytes.end - bytes.start;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("read"),
            size: size.max(4),
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        if size > 0 {
            encoder.copy_buffer_to_buffer(buffer, bytes.start, &staging, 0, size);
        }
        staging
    }

    pub async fn take<T: bytemuck::Pod>(&self, reading: Reading) -> Result<Vec<T>, String> {
        let slot = Arc::new(Mutex::new(Shot { value: None, waker: None }));
        let sender = slot.clone();
        let slice = reading.staging.slice(..);
        slice.map_async(wgpu::MapMode::Read, move |r| {
            if let Ok(mut s) = sender.lock() {
                s.value = Some(r);
                if let Some(w) = s.waker.take() {
                    w.wake();
                }
            }
        });
        #[cfg(not(target_arch = "wasm32"))]
        self.wait(Some(reading.index))?;
        #[cfg(target_arch = "wasm32")]
        let _ = reading.index;
        Oneshot(slot).await.map_err(|e| e.to_string())?;
        let data = slice.get_mapped_range();
        let out = bytemuck::cast_slice(&data).to_vec();
        drop(data);
        reading.staging.unmap();
        Ok(out)
    }

    pub async fn read_async<T: bytemuck::Pod>(
        &self,
        buffer: &wgpu::Buffer,
        elements: Range<usize>,
    ) -> Result<Vec<T>, String> {
        let size = std::mem::size_of::<T>() as u64;
        let bytes = elements.start as u64 * size..elements.end as u64 * size;
        let encoder = self.device.create_command_encoder(&Default::default());
        self.take(self.fetch(encoder, buffer, bytes)).await
    }

    pub fn length<T>(buffer: &wgpu::Buffer) -> usize {
        (buffer.size() / std::mem::size_of::<T>() as u64) as usize
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn take_now<T: bytemuck::Pod>(&self, reading: Reading) -> Result<Vec<T>, String> {
        pollster::block_on(self.take(reading))
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn read<T: bytemuck::Pod>(
        &self,
        buffer: &wgpu::Buffer,
        elements: Range<usize>,
    ) -> Result<Vec<T>, String> {
        pollster::block_on(self.read_async(buffer, elements))
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn read_all<T: bytemuck::Pod>(&self, buffer: &wgpu::Buffer) -> Result<Vec<T>, String> {
        self.read(buffer, 0..Gpu::length::<T>(buffer))
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn wait(&self, index: Option<wgpu::SubmissionIndex>) -> Result<(), String> {
        self.device
            .poll(wgpu::PollType::Wait { submission_index: index, timeout: None })
            .map(|_| ())
            .map_err(|e| e.to_string())
    }
}

struct Shot<T> {
    value: Option<T>,
    waker: Option<Waker>,
}

struct Oneshot<T>(Arc<Mutex<Shot<T>>>);

impl<T> Future for Oneshot<T> {
    type Output = T;

    fn poll(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<T> {
        let Ok(mut s) = self.0.lock() else {
            return Poll::Pending;
        };
        match s.value.take() {
            Some(v) => Poll::Ready(v),
            None => {
                s.waker = Some(cx.waker().clone());
                Poll::Pending
            }
        }
    }
}
