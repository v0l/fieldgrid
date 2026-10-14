# sims

Physics solvers that do not depend on any one app. They came out of
[gcad](https://github.com/v0l/gcad), [agentee](https://github.com/v0l/agentee) and
[antenna-toolbox](https://github.com/v0l/antenna-toolbox), which now build on them. Each solver
is its own crate, so an app pulls in only what it uses.

| crate | what it solves |
|---|---|
| `sims-field` | Rectilinear grids, triangle BVH, mesh rasterising, a wgpu compute harness, and a Yee-grid FDTD engine with CPML, lumped ports, Debye loss and near-to-far-field |
| `sims-flow` | A GPU lattice Boltzmann solver for air: D3Q19 with LES, coupled D3Q7 heat, turning rotors, drag and surface pressure on voxelised parts |
| `sims-fem` | Quadratic tetrahedral finite elements: elastic statics with contact, friction and bonded parts, heat conduction, modal (LOBPCG) and buckling analysis |
| `sims-tet` | Delaunay tetrahedral meshing of a closed surface mesh: boundary recovery, quality refinement, sliver removal, flips and smoothing |
| `sims-magnetic` | 2D magnetostatics for motors and actuators: a slice mesher with a turning airgap, nonlinear A_z with B-H curves, torque and flux linkage |
| `sims-circuit` | A SPICE-style circuit simulator: operating point, AC sweep, transient and noise, with diodes, BJTs, MOSFETs, op-amps and LDOs |
| `sims-quasistatic` | 2D transmission-line cross-sections (impedance, coupling, delay), conductor and dielectric loss with roughness and ENIG, and a 3D DC conduction grid, on the GPU or the CPU |
| `sims-antenna` | Antennas: method of moments for wires and surfaces, physical optics, Sommerfeld ground, characteristic modes, near fields, NEC decks, FDTD models from boxes and STL meshes |
| `sims-propagation` | Radio paths over terrain: ITM (Longley-Rice), ITU-R P.452, P.528 and P.2108, elevation tiles, path profiles and coverage maps |
| `sims-hf` | HF sky-wave circuits by ITU-R P.533, with P.372 radio noise |
| `sims-rf` | S-parameter networks: Touchstone, joining ports, noise, stability, step and TDR responses |
| `sims-progress` | Named progress stages that long solves report to one listener |

Pin a revision from git. A crate is found by name:

```toml
sims-fem = { git = "https://github.com/v0l/sims", rev = "..." }
sims-field = { git = "https://github.com/v0l/sims", rev = "..." }
```

## Using the mesher and FEM

`sims-tet` takes a closed triangle mesh as a `SurfaceMesh<S>`, where `S` says what each tagged
face really is. The mesher and the quadratic elements snap new boundary nodes onto it through the
`Surface` trait. gcad passes exact B-rep surfaces; `()` keeps the triangles as they are.

```rust
let tets = sims_tet::fill(surface, 2.5)?;
let mesh = sims_fem::Quadratic::from(&tets);
```

Then describe the supports and loads by face tag in an `elastic::Problem` and call
`elastic::solve`. `crates/sims-fem/tests/bar.rs` does this for a bar pulled at one end and
checks the stretch against FL/AE.

## GPUs

GPU solvers take an `Option<&sims_field::Gpu>` and fall back to the CPU without one. Open a
device once with `sims_field::Gpu::open_blocking` and pass it to each. In the browser the
blocking readbacks are not available, so `sims-quasistatic` runs on the CPU there.

## Progress

Every solver reports what it is doing through `sims-progress`: FDTD time steps, circuit
analyses, meshing, factoring and iterating, lattice steps, Newton on the magnetic field. Call
`sims_progress::listen` with a closure to receive the stages, or `listen_here` to receive
only those opened on the current thread when several solves run at once. With no listener
they do nothing. A new solver has to report too, with a test that listens for its stages.

## Tests

`cargo test --workspace` runs everything. Tests that need a GPU skip when no adapter opens; with
`--features sims-field/gpu-tests` the `sims-field` ones fail instead. CI runs them on Mesa's
software Vulkan driver:

```
VK_ICD_FILENAMES=/usr/share/vulkan/icd.d/lvp_icd.json cargo test --release --workspace --features sims-field/gpu-tests
```

The `sims-flow` validation cases (Poiseuille, the lid-driven cavity against Ghia, Schaefer-Turek,
heat conduction, turning walls) need a real GPU and a few minutes:
`cargo test --release -p sims-flow --features gpu-tests`.

What each solver is checked against is in `docs/validation.md`; open work is in `todo.md`.

## Licence

MIT
