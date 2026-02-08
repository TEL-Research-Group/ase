# ASE Psi4 Tutorial Map

Psi4-specific tutorials are not currently present in `doc/tutorials/`. Use these files as workflow scaffolding and adapt calculator sections to Psi4.

## Relevant tutorial sources
- `doc/tutorials/tutorials.rst`
  - States that most examples use one calculator but can be adapted to others.
- `doc/tutorials/tut03_vibrations/vibrations.rst`
  - Useful pattern for vibrational workflows; adapt calculator assignment to Psi4 where applicable.

## How to adapt
- Keep geometry/build/analysis flow from tutorial.
- Replace calculator construction with `ase.calculators.psi4.Psi4(...)`.
- Revalidate non-periodic assumptions and Psi4-supported keywords.

## Quick extraction commands
- Command: `rg -n -i "calculator|plugged|tutorial" doc/tutorials/tutorials.rst`
- Command: `rg -n -i "vibration|hessian|forces|calculator" doc/tutorials/tut03_vibrations/vibrations.rst`
