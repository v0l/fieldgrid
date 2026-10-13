# fieldgrid

Grids, mesh rasterising, a wgpu compute harness and an FDTD engine, shared by the field
solvers in [gcad](https://github.com/v0l/gcad), [agentee](https://github.com/v0l/agentee) and
[antenna-toolbox](https://github.com/v0l/antenna-toolbox). It runs natively and in the browser
through WebGPU.

It is not on crates.io. Pin a revision:

```toml
fieldgrid = { git = "https://github.com/v0l/fieldgrid", rev = "..." }
```

## Modules

| module | what it holds |
|---|---|
| `grid` | Rectilinear grids: uniform and graded mesh lines, smoothing, padding, nearest line and cell lookups |
| `bvh` | A triangle BVH with ray casts |
| `raster` | Meshes onto grids: `fill` marks samples inside a mesh and walls thinner than a cell, `solid_cells` and `solid_edges` do the same for closed solids, `sheet_edges` snaps zero-thickness sheets onto Yee edges, `closed` tells the two kinds apart |
| `gpu` | Opens a wgpu device (`Options` picks the backends, whether a software adapter is allowed and the storage buffer count), makes buffers and pipelines, and reads back with or without blocking |
| `fdtd` | A Yee-grid FDTD engine (below) |

## FDTD

`fdtd::Sim` holds a graded grid with CPML, per-cell materials with Debye loss, metal and
surface-impedance sheets on edges, lumped ports and RLC elements, plane field DFTs and a
near-to-far-field box. It steps on the GPU with `run` (async, so it works on WebGPU) or
`run_blocking`, or on the CPU with `run_cpu`. The CPU stepper is a port of every kernel and the
tests hold the two to float rounding. `ntff::Far` turns the recorded box into far fields,
radiation intensity and radiated power.

The engine knows nothing about boards or antennas. Each app builds the `Sim`: agentee from a
PCB stack-up, antenna-toolbox from boxes and STL meshes. `tests/fdtd.rs` builds small ones by
hand.

## Example

`examples/cube.rs` finds the cells inside a closed mesh:

```rust
let grid = Grid::uniform([0.0; 3], [1.0; 3], [20; 3], 0);
let inside = solid_cells(&Bvh::new(cube), &grid);
```

```
cargo run --example cube
```

## Tests

`cargo test` runs everything. Tests that need a GPU skip when no adapter opens; with
`--features gpu-tests` they fail instead. CI runs them on Mesa's software Vulkan driver:

```
VK_ICD_FILENAMES=/usr/share/vulkan/icd.d/lvp_icd.json cargo test --release --features gpu-tests
```

Parallel device creation can crash some Vulkan validation layers, so the FDTD tests take a
shared lock.

## Licence

MIT
