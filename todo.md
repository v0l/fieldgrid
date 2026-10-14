# todo

Open work on the solvers. Items came from the apps' own lists (`gcad/TODO.md`, `gcad/docs/flow.md`,
`gcad/docs/magnetic.md`, `gcad/docs/simulation.md`, `agentee/TODO.md`,
`antenna-toolbox/docs/validation.md`); each says where. Those lists still hold app-side work
(spec lines, file formats, CAD import) and their authors close their own items there.

## Progress (hard requirement)

Every solver entry point that can run longer than a moment reports through `sims-progress`:
`stage` for a whole run, `stage_after` for inner loops that are usually quick. Each crate has
a test that listens and checks for its stages (`a_run_reports_its_time_steps`,
`a_run_reports_its_analyses`, `meshing_and_solving_report_progress`); tests of delayed stages
listen with `listen_here_to_all`, which hears them from their start. New solvers follow the
same rule.

- A listener is global (`listen`) or for one thread (`listen_here`). A stage opened on a
  rayon worker reaches only the global one, so a thread-scoped listener (agentee's viewer)
  misses `field::solve` in `sims-magnetic` when it runs inside gcad's parallel pose sweep.

## sims-flow

- Smagorinsky puts eddy viscosity in the cells touching a wall (skin friction near 0.04
  against 0.003), so streamlined blades read badly: the Caradonna-Tung rotor
  (`gcad/examples/flow/caradonna/hover.gsim`, measured 507 N) reads 309 N and about 3x the
  torque at 13 mm cells, and its wake opens instead of contracting to 0.78 R. Replace it with
  WALE or Vreman in the shader in `crates/sims-flow/src/lattice.rs`, or add a wall function,
  then recheck `gcad/examples/flow/cube.gsim`, the heatsink and `square_cylinder_re22000`.
  Scaling it by 0.1 at the wall moved 20 mm cells from 110 N to 173 N; 0 diverges.
  (gcad TODO)
- Local grid refinement. The rotor gets 15 cells on the chord and under 2 through the blade
  at 13 mm, and thrust went 110 N at 20 mm to 309 N at 13 mm without settling. 5 mm over the
  rotator head is near the card's memory, and round members need 10 cells across.
  (gcad TODO, `docs/flow.md`)
- Fan benchmarks: the EAA FAN-01 axial fan in a short duct (geometry and LDA data at
  zenodo.org/records/10787093), and the fan wake against an OpenFOAM sliding-mesh transient
  instead of the frozen rotor in `gcad/tools/foam/fan.sh`. `gcad/examples/flow/fan.gsim`
  torque and thrust are unverified until then. (gcad TODO)
- A spinning part cannot carry heat: `spin` with `heat` bails in `gcad/src/sim/wind.rs`
  because the per-cell conductivity table of the thermal lattice would have to move with the
  rotor cells (`crates/sims-flow/src/rotor.rs`). (gcad TODO)
- Steady finite-volume RANS on the same grid, as a second answer. (`docs/flow.md`)
- f16 throughput is 79% of FluidX3D's; `precision half` is unchecked on the rotator.
  (`docs/flow.md`)

## sims-fem

- Parts bond wherever they come within a tenth of their element size, so a cover with a
  small clearance is glued to the part it covers: in `~/git/rotator/sims/rotator.gsim`
  `cover_el` bonds to `wheel_el` at 29,830 nodes and `cover_az` to `wheel_az` at 5,049, and
  the solve goes from 592 to 1,372 iterations (that the clearance is what gets bridged is
  unverified; measure the gap in `parts/gearcover.gcad`). Bond only faces that touch, or take
  a list of bonded pairs. The distance is set in `gcad/src/sim/statics.rs`, the bonding is
  `merge` in `crates/sims-fem/src/lib.rs`. (gcad TODO)
- The static solve runs on the CPU only; a GPU solve waits on a multigrid preconditioner.
  (`gcad/docs/simulation.md`, phase 3)

## sims-magnetic

- Not built: `phase` and `drive` (sinusoidal phase currents that follow the rotor through a
  sweep) and a `kt` reading; force on a sliding part (weighted Maxwell stress); eddy and core
  losses (2D time-harmonic A_z, Steinmetz/Bertotti) feeding heat; 3D magnetics on the tet mesh
  checked against TEAM problem 13. (gcad TODO, `docs/magnetic.md`)

## sims-circuit

- No Pst flicker readout (EN 61000-3-3): it needs the IEC 61000-4-15 flickermeter on the
  supply voltage over minutes of transient. (agentee TODO)

## sims-hf and sims-propagation

- P.533 (`sims-hf`): the mirror-reflection height term G omits the `90.47·xr` term the Recommendation has
  for `xr <= 3.7` (ITURHFProp omits it too, and the port keeps it out so the two agree), and only
  analogue circuit reliability (BCR) is reported. P.452 (`sims-propagation`) takes ΔN and N₀ by hand instead of from
  the ITU digital maps, which may not be redistributed. (`antenna-toolbox/docs/validation.md`)

## Not in sims yet

- `sims-rf` has its own `Cx` complex type; everything else uses `num-complex`, and
  antenna-toolbox `crates/rf` keeps its own network maths.
