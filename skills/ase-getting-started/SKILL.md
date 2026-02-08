---
name: ase-getting-started
description: This skill should be used when users ask about getting started in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Getting Started

## High-Signal Playbook

### Route the request
- Use this skill for first-run setup, backend readiness checks, and first successful energy/force evaluation.
- Route detailed backend parameterization to `ase-calculators`.
- Route Psi4-specific setup and behavior questions to `ase-psi4-interface`.
- Route workflow questions (MD/NEB/QM-MM loops) to `ase-simulation-workflows`.
- Route scripting/database/constraint details to `ase-api-and-scripting`.

### Triage questions
- Which backend calculator is being used (for example VASP, SIESTA, OpenMX, ONETEP, KIM, DFTD3)? (`doc/ase/calculators/vasp.rst`, `doc/ase/calculators/siesta.rst`, `doc/ase/calculators/openmx.rst`, `doc/ase/calculators/onetep.rst`, `doc/ase/calculators/kim.rst`, `doc/ase/calculators/dftd3.rst`)
- Is the runtime command configured (environment variable or calculator `command`)? (`doc/ase/calculators/vasp.rst`, `doc/ase/calculators/siesta.rst`, `doc/ase/calculators/openmx.rst`, `doc/ase/calculators/onetep.rst`)
- Are pseudopotentials/model files discoverable from configured paths? (`doc/ase/calculators/vasp.rst`, `doc/ase/calculators/siesta.rst`, `doc/ase/calculators/openmx.rst`)
- Is this a one-shot single-point test or repeated ASE-driven steps where socket mode is better? (`doc/ase/calculators/socketio/socketio.rst`)
- Is GUI inspection needed for quick sanity checks? (`doc/ase/gui/basics.rst`)

### Canonical workflow
1. Confirm calculator availability and local configuration (use `ase info --calculators` when relevant to socket/backends) (`doc/ase/calculators/socketio/socketio.rst`).
2. Set backend command/pseudopotential variables for the chosen code (for example `ASE_VASP_COMMAND`, `ASE_SIESTA_COMMAND`, `OPENMX_DFT_DATA_PATH`, `ASE_ONETEP_COMMAND`) (`doc/ase/calculators/vasp.rst`, `doc/ase/calculators/siesta.rst`, `doc/ase/calculators/openmx.rst`, `doc/ase/calculators/onetep.rst`).
3. Run a minimal structure with one calculator and call `get_potential_energy()` to verify end-to-end execution (`doc/ase/calculators/kim.rst`).
4. If running iterative ASE loops with external DFT engines, consider Socket I/O client/server mode to reduce per-step startup overhead (`doc/ase/calculators/socketio/socketio.rst`).
5. Inspect geometry/path quickly with `ase gui` before scaling up (`doc/ase/gui/basics.rst`).

### Minimal working examples
```python
from ase.lattice.cubic import FaceCenteredCubic
from ase.calculators.kim.kim import KIM

atoms = FaceCenteredCubic(symbol='Ar', latticeconstant=5.25, size=(1, 1, 1))
atoms.calc = KIM("ex_model_Ar_P_Morse_07C")
print(atoms.get_potential_energy())
```

```bash
ase gui x.traj@-10:
ase gui -n -1 *.traj
```

### Pitfalls
- Missing backend command variables is the most common startup failure for file-IO calculators (`doc/ase/calculators/vasp.rst`, `doc/ase/calculators/siesta.rst`, `doc/ase/calculators/openmx.rst`, `doc/ase/calculators/onetep.rst`).
- DFT-D3 alone is a dispersion correction, not a complete dynamics-ready potential (`doc/ase/calculators/dftd3.rst`).
- DFT-D3 treats 1D/2D PBC systems as fully 3D periodic (`doc/ase/calculators/dftd3.rst`).
- Socket sessions should be closed cleanly (`with` context or `calc.close()`) to avoid hanging connections (`doc/ase/calculators/socketio/socketio.rst`).
- Gromacs wrapper is documented as slow and mainly useful for QM/MM/testing rather than pure MM production (`doc/ase/calculators/gromacs.rst`).

### Convergence/validation checklist
- Confirm `atoms.get_potential_energy()` and (if needed) `atoms.get_forces()` return finite values.
- Verify external command and pseudo/model paths resolve to real files.
- Confirm the intended backend appears in `ase info --calculators`.
- For socket workflows, run one short loop and verify graceful termination.

## Scope
- Handle questions about initial setup, quickstarts, and core concepts.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/calculators/harmonic.rst`
- `doc/ase/calculators/octopus.rst`
- `doc/ase/calculators/dftb.rst`
- `doc/ase/gui/basics.rst`
- `doc/ase/calculators/vasp.rst`
- `doc/ase/calculators/siesta.rst`
- `doc/ase/calculators/crystal.rst`
- `doc/ase/calculators/openmx.rst`
- `doc/ase/calculators/plumed.rst`
- `doc/ase/calculators/onetep.rst`
- `doc/ase/calculators/gromacs.rst`
- `doc/ase/calculators/FHI-aims.rst`

## Workflow
- Start with the primary references above.
- If details are missing, inspect `references/doc_map.md` for the complete topic document list.
- Use tutorials/examples as executable usage patterns when available.
- Use tests as behavior or regression references when available.
- If ambiguity remains after docs, inspect `references/source_map.md` and start with the ranked source entry points.
- Cite exact documentation file paths in responses.

## Tutorials and examples
- `doc/tutorials`

## Test references
- `ase/test`

## Optional deeper inspection
- `ase`

## Source entry points for unresolved issues
- `ase/calculators/dftd3.py`
- `ase/calculators/dftb.py`
- `ase/calculators/siesta/siesta_lrtddft.py`
- `ase/calculators/vasp/__init__.py`
- `ase/calculators/siesta/__init__.py`
- `ase/calculators/openmx/__init__.py`
- `ase/calculators/kim/__init__.py`
- `ase/calculators/exciting/__init__.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
