# Validation

What each solver has been checked against. Tests named with a path are in this repo; the rest
are in the app the solver came from, which drives it through its own files.

## sims-fem and sims-tet

Against gmsh and scikit-fem on the same STEP exports (`gcad/tools/xcheck/run.sh`), sharing no code:

| case | sims | reference |
|---|---|---|
| bracket, mesh 1.2: peak displacement | 0.003337 mm | 0.003326 mm |
| same, `top` probe | 0.001786 mm, 3.762 MPa | 0.001783 mm, 3.788 MPa |
| bracket fixed at its base, first four frequencies | 43804, 43854, 69358, 99613 Hz | 43911, 43979, 69490, 99741 Hz |
| heatsink: hottest, coolest | 64.833, 63.935 °C | 64.832, 63.935 °C |
| rod in tube between plates, rod cooled 100 °C | 142.2, 28.12 MPa | 143.3, 27.56 MPa |
| 10 x 10 x 200 steel column, buckling factors | 10.285, 10.285, 91.08 | 10.285, 10.285, 91.09 |
| 20 x 20 x 3 angle, 150 long, buckling factors | 36.48, 102.6, 211.7 | 36.48, 102.5, 211.4 |

`crates/sims-fem/tests/bar.rs`: a bar meshed from plain triangles stretches within 1% of FL/AE.
Contact and friction rest on gcad's stack tests.

## sims-flow

| case | sims | reference |
|---|---|---|
| Poiseuille between plates, round pipe | within 0.2%, 0.15% | exact |
| lid-driven cavity, Re 100 and 1000 | within 0.004 and 0.015 of the lid speed | Ghia et al. |
| Schaefer-Turek cylinder, Re 20, 60 cells | drag 5.602, pressure drop 2.9429 | 5.5795, 2.9375 |
| Schaefer-Turek cylinder, Re 100 | drag 3.229, Strouhal 0.302 | 3.22-3.24, 0.295-0.305 |
| square cylinder, Re 22,000 | drag 2.20, Strouhal 0.145 | Lyn 2.1, 0.132 (LES spread 1.66-2.77) |
| channel between heated walls, Re 50 | Nu 8.232, 8.239 | 8.235 |
| sphere, Re 50, 20 / 40 cells (gcad) | 1.6265 / 1.6150 | OpenFOAM, same box, 1.606 |
| cube, 10 m/s, 20 cells (gcad) | Cd 1.01 | Hoerner 1.05 |
| heatsink, 0.5 / 0.33 / 0.25 mm (gcad) | 1.925 / 1.854 / 1.845 K/W | OpenFOAM 1.914, Teertstra 1.887 |

The first six are `crates/sims-flow/tests/flow.rs`, run with `--features gpu-tests`. The
comparison with FluidX3D and OpenFOAM on the same boxes is in `gcad/docs/flow.md`.

## sims-magnetic

`crates/sims-magnetic/tests/magnetic.rs`: a round wire's field within 0.9% RMS of the surface
field; a round magnet holds 0.599 T of 0.6 inside and its dipole field outside within 0.7%; a
saturating M-19 toroid carries the flux of its B-H curve within 1%; the torque a wire puts on a
magnet within 1% at every angle; the torque of a 12 slot 10 pole motor against
xfemm at seven rotor angles.

## sims-antenna, sims-propagation, sims-hf

| model | reference | test | agreement |
|---|---|---|---|
| Longley-Rice ITM | NTIA ITM v1.4, SPLAT! 1.4.2 | `sims_propagation::itm` | 0.01 dB on NTIA p2p cases |
| ITU-R P.528 | NTIA p528 | `sims_propagation::p528` | 0.01 dB on 400 paths |
| ITU-R P.452-18 | Py452 validation set | `sims_propagation::p452` | 2e-7 dB on 595 cases |
| ITU-R P.2108 | NTIA p2108 test data | `sims_propagation::p2108` | every row |
| ITU-R P.533-14, P.372 | ITURHFProp 14.2 | `crates/sims-hf/tests/p533.rs` | 0.005 on 1200 cases |
| near fields | nec2c NE/NH cards | `crates/sims-antenna/tests/validate.rs` | 1% |
| characteristic modes | scipy `eig(X, R)` | `crates/sims-antenna/tests/validate.rs` | 1e-3 |
| FDTD dipole | energy balance, MoM | `crates/sims-antenna/tests/fdtd.rs` | 1% power, 2.15 dBi |
| FDTD thin-wire dipole, radius 0.06 to 0.45 mm | the moment method, same radius | `crates/sims-antenna/tests/thin_wire.rs` | within 1.2% at 1.5 mm cells, 2.3% at 3 mm |
| FDTD microstrip patch, boxes and meshes | Sheen et al. 1990, measured 7.5 GHz | `crates/sims-antenna/tests/fdtd.rs` | within 1.5% |

## sims-field

`crates/sims-field/tests/fdtd.rs`: an empty cavity rings at its analytic modes; every kernel
steps the same on the CPU and the GPU (4e-7; WebGPU 2.7e-4); thin wires at both ends of their
radius range ring in a closed cavity for 40,000 steps without growing.

Thin wires scale the permeability of the four magnetic components circling the wire by 1/f and
the permittivity of the electric edges leaving it by f, with f = 2 / ln(cell / radius), so the
first ring of cells carries the inductance and capacitance of the wire it stands for. Scaling
both sides keeps the update symmetric; scaling only the magnetic side (Taflove's contour form)
grew without bound after about 8,000 steps. The time step shrinks by the square root of the
largest f or 1/f, which limits the radius to between a 55th and half of a cell.

## sims-circuit

`crates/sims-circuit/tests`: a divider, an RC corner (3.010 dB, 45°) and charge (63%); a
1 V/us follower ramps at its slew rate; harmonic currents against EN 61000-3-2 on a known third
harmonic; a line's level against CISPR 11; the 1.2/50 us and 8/20 us surges against the
standard's front and duration; 5/50 ns burst pulses; a varistor's power law and its clamp on a
2 ohm surge. agentee's `crates/agentee-sim/tests/circuit.rs` adds circuits from schematics: diodes,
BJTs, MOSFETs, op-amps at their rails, LDOs, noise, a LISN against the hand calculation.

## sims-quasistatic

Stripline against Cohn's exact solution within 0.5%, edge-coupled pairs within 1%, microstrip
against Hammerstad-Jensen within 0.5% (in agentee), conductor loss against Wheeler's formula in
Pozar within 3%, the TEM dielectric loss of a homogeneous line within 0.2%, and the roughness and
ENIG surface impedance against their published limits.

## sims-logic

`crates/sims-logic/src/engine.rs`: every primitive against its truth table, unknowns that
propagate only when they decide the output, strong drivers over weak and contention, setup,
hold, recovery and removal against the clock, and VCD buses grouped as the viewer shows them.
agentee's `crates/agentee-sim/tests/logic.rs` runs its counter and I2C examples from schematics.

## sims-frame

`crates/sims-frame/tests/frame.rs`: a cantilever's tip deflection both ways and its twist, a
simply supported beam under a uniform load and the same beam made of fixed members with pinned
ends, all to 1e-9 of the beam formulas; a fixed portal's sway within 0.2% of slope deflection
(the rest is the columns' axial give); Euler's load for a pinned and a cantilevered column and
a cantilever's first two bending frequencies within 1e-4 on ten and twenty elements. A
three-storey, two-by-two-bay concrete frame with floor masses matches OpenSees
(`elasticBeamColumn`, `tools/opensees_frame.py`) within 1e-6 on displacements, a base moment
and the first six periods.
