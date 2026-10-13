# fieldgrid

Shared pieces of the field solvers in gcad, agentee and antenna-toolbox:

- `grid`: rectilinear grids and graded mesh lines
- `bvh` and `raster`: triangle meshes and filling grid samples inside them
- `gpu`: a wgpu compute harness that opens a device, makes buffers and pipelines, and reads back on native and in the browser

`cargo test --features gpu-tests` fails instead of skipping when no adapter opens.
