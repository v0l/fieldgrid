# fieldgrid

Shared pieces of the field solvers in gcad, agentee and antenna-toolbox:

- `grid`: rectilinear grids and graded mesh lines
- `bvh` and `raster`: triangle meshes and filling grid samples inside them
- `gpu`: a wgpu compute harness that opens a device, makes buffers and pipelines, and reads back on native and in the browser
- `fdtd`: a Yee-grid FDTD engine on graded grids with CPML, lumped ports and elements, Debye dielectric loss, surface-impedance sheets, plane field DFTs and near-to-far-field patches. It steps on the GPU (`run`, async, native and WebGPU) or the CPU (`run_cpu`), and the two agree to float rounding

`cargo test --features gpu-tests` fails instead of skipping when no adapter opens.
