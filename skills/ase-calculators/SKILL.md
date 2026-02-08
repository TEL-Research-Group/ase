---
name: ase-calculators
description: This skill should be used when users ask about calculators in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Calculators

## High-Signal Playbook

### Route the request
- Use this skill for selecting/configuring calculators and composing wrapper calculators.
- Route Psi4-specific integration and troubleshooting questions to `ase-psi4-interface`.
- Route structure/model-definition questions to `ase-inputs-and-modeling`.
- Route runtime loop orchestration (MD/NEB/QM-MM runs) to `ase-simulation-workflows`.

### Triage questions
- Is the target a pure ASE calculator, external-code wrapper, or composite wrapper (`Checkpoint`, `Mixing`, `SocketIO`, `QMMM`)? (`doc/ase/calculators/calculators.rst`, `doc/ase/calculators/checkpointing.rst`, `doc/ase/calculators/qmmm.rst`)
- Which external executables and pseudopotentials are required?
- Is restart/rollback behavior needed?
- Is the workload dominated by expensive SCF steps or cheap classical steps?

### Canonical workflow
1. Select calculator class and verify prerequisites (binary, config, pseudo/model files) (`doc/ase/calculators/calculators.rst`).
2. Attach calculator and run a single-point smoke test (`doc/ase/calculators/calculators.rst`).
3. Add wrapper calculators only when the use case justifies them (checkpointing, mixing, socket, QM/MM) (`doc/ase/calculators/calculators.rst`, `doc/ase/calculators/checkpointing.rst`, `doc/ase/calculators/qmmm.rst`).
4. For long expensive workflows, prefer checkpointing/robust restart boundaries (`doc/ase/calculators/checkpointing.rst`).
5. Validate energy/forces/stress consistency before production sweeps.

### Minimal working examples
```python
from ase.io import read
from ase.calculators.abinit import Abinit

atoms = read('molecule.xyz')
atoms.calc = Abinit(...)
print(atoms.get_potential_energy())
```

```python
from ase.calculators.checkpoint import CheckpointCalculator

cp_calc = CheckpointCalculator(calc)
atoms.calc = cp_calc
print(atoms.get_potential_energy())
```

### Pitfalls
- Calling energy/forces before setting `atoms.calc` raises a runtime error (`doc/ase/calculators/calculators.rst`).
- DFT-D3 corrections are not standalone dynamics potentials (`doc/ase/calculators/dftd3.rst`).
- `CheckpointCalculator` can create huge databases for classical MD and is not recommended there (`doc/ase/calculators/checkpointing.rst`).
- Command construction should include launcher/binary but avoid unnecessary extra flags unless intentional (`doc/ase/calculators/calculators.rst`).
- Composite workflows (QM/MM, mixing) require explicit validation of coupling assumptions (`doc/ase/calculators/qmmm.rst`).

### Convergence/validation checklist
- Confirm `energy`, `forces`, and (if needed) `stress` are available and finite.
- Verify restart path by re-running from saved state.
- Benchmark wrapper overhead versus plain calculator on a small test case.
- Ensure calculator configuration is captured in script/config for reproducibility.

## Scope
- Handle questions about documentation grouped under the 'calculators' theme.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/calculators/others.rst`
- `doc/ase/calculators/qmmm.rst`
- `doc/ase/calculators/checkpointing.rst`
- `doc/ase/calculators/calculators.rst`
- `doc/ase/calculators/ace.rst`
- `doc/ase/calculators/abinit.rst`
- `doc/ase/calculators/test.rst`
- `doc/ase/calculators/mopac.rst`
- `doc/ase/calculators/mixing.rst`
- `doc/ase/calculators/loggingcalc.rst`
- `doc/ase/calculators/lammpslib.rst`
- `doc/ase/calculators/lammps.rst`

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
- `ase/calculators/lammpslib.py`
- `ase/calculators/lammps/__init__.py`
- `ase/calculators/demon/__init__.py`
- `ase/calculators/qmmm.py`
- `ase/calculators/mopac.py`
- `ase/calculators/mixing.py`
- `ase/calculators/loggingcalc.py`
- `ase/calculators/fleur.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
