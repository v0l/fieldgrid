# sims

Physics solvers that do not depend on any one app. They came out of
[gcad](https://github.com/v0l/gcad), [agentee](https://github.com/v0l/agentee) and
[antenna-toolbox](https://github.com/v0l/antenna-toolbox), which now build on them. Each solver
is its own crate, so an app pulls in only what it uses.

| crate | what it solves |
|---|---|
| `fieldgrid` | Rectilinear grids, triangle BVH, mesh rasterising, a wgpu compute harness, and a Yee-grid FDTD engine with CPML, lumped ports, Debye loss and near-to-far-field |
| `sims-fem` | Quadratic tetrahedral finite elements: elastic statics with contact, friction and bonded parts, heat conduction, modal (LOBPCG) and buckling analysis |
| `sims-tet` | Delaunay tetrahedral meshing of a closed surface mesh: boundary recovery, quality refinement, sliver removal, flips and smoothing |
| `sims-circuit` | A SPICE-style circuit simulator: operating point, AC sweep, transient and noise, with diodes, BJTs, MOSFETs, op-amps and LDOs |
| `sims-quasistatic` | 2D transmission-line cross-sections (impedance, coupling, delay), conductor and dielectric loss with roughness and ENIG, and a 3D DC conduction grid, on the GPU or the CPU |
| `sims-rf` | S-parameter networks: Touchstone, joining ports, noise, stability, step and TDR responses |
| `sims-progress` | Named progress stages that long solves report to one listener |

Pin a revision from git. A crate is found by name:

```toml
sims-fem = { git = "https://github.com/v0l/sims", rev = "..." }
fieldgrid = { git = "https://github.com/v0l/sims", rev = "..." }
```

## Using the mesher and FEM

`sims-tet` takes a closed triangle mesh as a `SurfaceMesh<S>`, where `S` says what each tagged
face really is. The mesher and the quadratic elements snap new boundary nodes onto it through the
`Surface` trait. gcad passes exact B-rep surfaces; `()` keeps the triangles as they are.

```rust
let tets = sims_tet::fill(surface, 2.5, &|_| {})?;
let mesh = sims_fem::Quadratic::from(&tets);
```

Then describe the supports and loads by face tag in an `elastic::Problem` and call
`elastic::solve`. `crates/sims-fem/tests/bar.rs` does this for a bar pulled at one end and
checks the stretch against FL/AE.

## GPUs

GPU solvers take an `Option<&fieldgrid::Gpu>` and fall back to the CPU without one. Open a
device once with `fieldgrid::Gpu::open_blocking` and pass it to each. In the browser the
blocking readbacks are not available, so `sims-quasistatic` runs on the CPU there.

## Progress

`sims-fem` reports stages such as factoring and iterating through `sims-progress`. Call
`sims_progress::listen` with a closure to receive them; with no listener they do nothing.

## Tests

`cargo test --workspace` runs everything. Tests that need a GPU skip when no adapter opens; with
`--all-features` the `fieldgrid` ones fail instead. CI runs them on Mesa's software Vulkan
driver:

```
VK_ICD_FILENAMES=/usr/share/vulkan/icd.d/lvp_icd.json cargo test --release --workspace --all-features
```

## Licence

MIT
