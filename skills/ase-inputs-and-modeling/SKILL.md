---
name: ase-inputs-and-modeling
description: This skill should be used when users ask about inputs and modeling in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Inputs and Modeling

## High-Signal Playbook

### Route the request
- Use this skill for structure definition, cell/PBC setup, model parameterization, and sampling choices.
- Route calculator-engine-specific flags and external-code setup to `ase-calculators`.
- Route runtime protocol questions (MD/NEB/QM-MM execution loops) to `ase-simulation-workflows`.

### Triage questions
- Is the system molecular, slab, wire, or bulk periodic crystal? (`doc/ase/atoms.rst`, `doc/ase/spacegroup/spacegroup.rst`, `doc/ase/cluster/cluster.rst`)
- Are cell vectors and PBC flags set consistently with intended physics? (`doc/ase/atoms.rst`)
- Which basis/force-field/model choices are required (for example Tersoff, Espresso pseudopotentials)? (`doc/ase/calculators/tersoff.rst`, `doc/ase/calculators/espresso.rst`)
- What Brillouin-zone sampling strategy is planned (Gamma-only, Monkhorst-Pack, path)? (`doc/ase/dft/kpoints.rst`)
- Do you need robust IO interoperability before running expensive jobs? (`doc/ase/io/io.rst`)

### Canonical workflow
1. Build the structure using `Atoms`, `ase.spacegroup.crystal`, or cluster builders as appropriate (`doc/ase/atoms.rst`, `doc/ase/spacegroup/spacegroup.rst`, `doc/ase/cluster/cluster.rst`).
2. Set cell and PBC explicitly, and decide whether coordinate scaling is required when changing the cell (`doc/ase/atoms.rst`).
3. Select calculator/model parameters and attach them to the structure (`doc/ase/calculators/espresso.rst`, `doc/ase/calculators/tersoff.rst`).
4. Choose and test k-point sampling (Monkhorst-Pack or symmetry path) (`doc/ase/dft/kpoints.rst`).
5. Run a short geometry optimization and inspect forces/energy trend (`doc/ase/optimize.rst`).
6. Verify IO round-trips (write/read) for the chosen exchange format before large sweeps (`doc/ase/io/io.rst`).

### Minimal working examples
```python
from ase import Atoms
import numpy as np

atoms = Atoms('N3', [(0, 0, 0), (1, 0, 0), (0, 0, 1)])
atoms.set_cell(2 * np.identity(3))
atoms.set_pbc((True, True, False))
```

```python
from ase.dft.kpoints import monkhorst_pack
print(monkhorst_pack((4, 1, 1)))
```

### Pitfalls
- `set_cell()` does not move atoms unless `scale_atoms=True` is passed (`doc/ase/atoms.rst`).
- Inconsistent PBC flags versus geometry type leads to incorrect boundary physics (`doc/ase/atoms.rst`).
- Using unconverged k-point meshes can dominate total-energy and force errors (`doc/ase/dft/kpoints.rst`).
- Tersoff parameters must be complete and physically consistent for each interaction tuple (`doc/ase/calculators/tersoff.rst`).
- IO conversions can lose uncommon format features; verify assumptions explicitly (`doc/ase/io/io.rst`).

### Convergence/validation checklist
- Confirm geometry is physically plausible after any cell transform.
- Perform k-point sensitivity checks for periodic systems.
- Ensure optimization reaches target `fmax` with stable energy trend.
- Round-trip at least one representative structure through the chosen IO format.

## Scope
- Handle questions about inputs, system setup, models, and physical parameterization.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/atoms.rst`
- `doc/ase/optimize.rst`
- `doc/ase/lattice.rst`
- `doc/ase/dft/kpoints.rst`
- `doc/ase/calculators/tersoff.rst`
- `doc/tutorials/wannier/wannier.rst`
- `doc/ase/thermochemistry/thermochemistry.rst`
- `doc/ase/spacegroup/spacegroup.rst`
- `doc/ase/calculators/espresso.rst`
- `doc/ase/cluster/cluster.rst`
- `doc/ase/geometry.rst`
- `doc/ase/io/io.rst`

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
- `ase/dft/band_structure.py`
- `ase/dft/wannier.py`
- `ase/dft/kpoints.py`
- `ase/calculators/tersoff.py`
- `ase/calculators/espresso.py`
- `ase/calculators/demonnano.py`
- `ase/calculators/vasp/create_input.py`
- `ase/geometry/__init__.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
