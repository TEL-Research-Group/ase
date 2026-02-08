---
name: ase-psi4-interface
description: This skill should be used when users ask about using Psi4 through ASE (setup, parameters, limitations, threading/scratch behavior, energy/force calculations, reading psi4-calc outputs, or troubleshooting ASE+Psi4 workflows), with answers grounded in ASE documentation, source code, tests, and relevant tutorial patterns.
---

# ase: Psi4 Interface

## High-Signal Playbook

### Route the request
- Use this skill for Psi4-specific ASE questions.
- Route generic calculator-selection questions to `ase-calculators`.
- Route non-Psi4 workflow orchestration (MD/NEB/QM/MM) to `ase-simulation-workflows`.

### Triage questions
- Is the user asking about setup (installation/scratch/threads) or runtime behavior?
- Is the system non-periodic? (Psi4 calculator rejects periodic `kpts` usage.)
- Are they using `method`/`basis` correctly, and avoiding unsupported ASE-style kwargs (`xc`, `smearing`, `nbands`)?
- Do they need only energy, or both energy and forces?
- Are they trying to read an existing `*.dat` file that may not have ASE’s embedded metadata block?

### Canonical workflow
1. Confirm environment and installation assumptions from `references/doc_map.md`.
2. Build a molecule/non-periodic `Atoms` object and attach `Psi4(...)`.
3. Set `method`, `basis`, optional `memory`, optional `num_threads`, and scratch path (`PSI_SCRATCH`) if needed.
4. Run `get_potential_energy()` and optionally `get_forces()`.
5. If needed, use `calc.psi4` for advanced Psi4 API calls (for example frequency workflows).
6. For restart/inspect flows, read with `calc.read(label)` only when the `label.dat` file was produced by ASE Psi4 calculator (contains `!ASE Information` JSON block).

### Minimal working examples
```python
from ase.build import molecule
from ase.calculators.psi4 import Psi4

atoms = molecule('H2O')
atoms.calc = Psi4(method='b3lyp', basis='6-311g_d_p_', memory='500MB')
print(atoms.get_potential_energy())
print(atoms.get_forces())
```

```python
from ase.calculators.psi4 import Psi4

calc = Psi4()
calc.read('psi4-calc')  # reads psi4-calc.dat written by ASE Psi4 calculator
print(calc.results['energy'])
```

### Pitfalls
- Psi4 calculator is non-periodic; `kpts` raises `InputError`.
- `xc` is rejected; use `method` instead.
- `smearing` and `nbands` are rejected.
- `method='LDA'` is remapped internally to `svwn`.
- If charge is set without multiplicity, multiplicity defaults to `1` with a warning.
- `read()` fails for Psi4 output files that do not contain ASE’s `!ASE Information` block.
- `command`/shell-launch semantics are not the standard FileIO-calculator path here; ASE uses the Psi4 Python API.

### Convergence/validation checklist
- Validate method/basis with a small molecule first.
- Check energy/force consistency (optionally against numerical forces).
- Confirm thread count and scratch directory are intentional for the machine.
- Keep the generated `label.dat` file for reproducibility and ASE-readable restart metadata.

## Scope
- Handle ASE+Psi4 interface questions: setup, parameters, implementation behavior, output reading, and troubleshooting.
- Keep answers grounded in docs first, then source/tests when behavior details are needed.

## Primary documentation references
- `doc/ase/calculators/psi4.rst`
- `doc/ase/calculators/calculators.rst`

## Supporting references
- `references/doc_map.md`
- `references/tutorial_map.md`
- `references/source_map.md`
- `references/test_map.md`

## Workflow
- Start with documentation for user-facing setup and parameters.
- Use tutorials to shape practical workflow scaffolding.
- Resolve ambiguities with source-level behavior from `ase/calculators/psi4.py`.
- Use tests to confirm expected numerical/IO behavior and known limitations.
- Cite exact file paths when answering.

## Tutorials and examples
- `doc/tutorials`

## Test references
- `ase/test`

## Optional deeper inspection
- `ase`

## Source entry points for unresolved issues
- `ase/calculators/psi4.py`
- `ase/test/calculator/psi4/test_psi4_HF_3_21G.py`
- `ase/test/factories.py`
- `ase/test/calculator/test_command.py`
