---
name: ase-simulation-workflows
description: This skill should be used when users ask about simulation workflows in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Simulation Workflows

## High-Signal Playbook

### Route the request
- Use this skill for execution loops: MD, QM/MM runs, and NEB-style progression/analysis.
- Route static model-building questions to `ase-inputs-and-modeling`.
- Route backend installation/command specifics to `ase-calculators`.

### Triage questions
- Which workflow type is needed: MD, NEB pathway, or QM/MM? (`doc/ase/md.rst`, `doc/tutorials/dissociation.rst`, `doc/tutorials/qmmm/qmmm.rst`)
- What ensemble and timestep are planned for MD? (`doc/ase/md.rst`)
- Is restart/resume needed for long jobs? (`doc/tutorials/neb/diffusion.rst`)
- For QM/MM, what partition and coupling model are required? (`doc/ase/calculators/qmmm.rst`, `doc/tutorials/qmmm/qmmm.rst`)
- For LAMMPS integration, are current wrapper limits acceptable? (`doc/ase/calculators/lammpsrun.rst`)

### Canonical workflow
1. Start from a relaxed structure and define constraints/region partitioning when relevant (`doc/tutorials/constraints/diffusion.rst`, `doc/tutorials/qmmm/qmmm.rst`).
2. Pick the dynamics/path driver (`VelocityVerlet`, NEB, or QM/MM calculator) (`doc/ase/md.rst`, `doc/tutorials/dissociation.rst`, `doc/ase/calculators/qmmm.rst`).
3. Set timestep/thermostat choices and output intervals intentionally (`doc/ase/md.rst`).
4. Run a short pilot to validate stability, then continue production (`doc/ase/md.rst`).
5. Analyze trajectories/barriers with GUI and NEB tooling where applicable (`doc/tutorials/neb/diffusion.rst`).

### Minimal working examples
```python
from ase import units
from ase.md.verlet import VelocityVerlet

dyn = VelocityVerlet(atoms, dt=5.0 * units.fs, trajectory='md.traj', logfile='md.log')
dyn.run(1000)
```

```python
from ase.calculators.qmmm import EIQMMM, LJInteraction
from ase.calculators.tip3p import epsilon0, sigma0

lj = LJInteraction({'OO': (epsilon0, sigma0)})
atoms.calc = EIQMMM([0, 1, 2], QMCalculator(...), MMCalculator(...), interaction=lj)
```

### Pitfalls
- Too-large timestep causes instability and energy drift/blow-up (`doc/ase/md.rst`).
- MD APIs now prefer `temperature_K`; older `temperature` usage may warn (`doc/ase/md.rst`).
- Writing/logging every step can create unnecessary IO overhead (`doc/ase/md.rst`).
- QM/MM force-based schemes require explicit buffer-width testing (`doc/ase/calculators/qmmm.rst`).
- LAMMPSrun is documented as a thin wrapper; feature expectations should be calibrated (`doc/ase/calculators/lammpsrun.rst`).

### Convergence/validation checklist
- Check NVE energy drift on a short pilot run.
- Verify temperature/pressure behavior matches intended ensemble.
- Confirm trajectory/log intervals are sufficient for post-analysis without exploding storage.
- For NEB, verify barrier/profile consistency in GUI (`Tools -> NEB`).
- For QM/MM, check partition consistency and force continuity at boundaries.

## Scope
- Handle questions about simulation setup, execution flow, and runtime controls.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/xrdebye.rst`
- `doc/ase/md.rst`
- `doc/tutorials/qmmm/qmmm.rst`
- `doc/ase/calculators/lammpsrun.rst`

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
- `ase/calculators/lammpsrun.py`
- `ase/calculators/qmmm.py`
- `ase/io/lammpsrun.py`
- `ase/calculators/castep.py`
- `ase/calculators/exciting/runner.py`
- `ase/ga/pbs_queue_run.py`
- `ase/cli/run.py`
- `ase/calculators/__init__.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
