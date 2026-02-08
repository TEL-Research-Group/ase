---
name: ase-parallel-hpc
description: This skill should be used when users ask about parallel and hpc in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Parallel and HPC

## High-Signal Playbook

### Route the request
- Use this skill for MPI/threading strategy, per-rank IO, and image-parallel NEB patterns.
- Route calculator-specific parallel flags and backend setup details to `ase-calculators`.
- Route Psi4-specific setup/API constraints to `ase-psi4-interface`.
- Route path setup/analysis details to `ase-simulation-workflows`.

### Triage questions
- Is parallelism over atoms, over images, or over independent jobs? (`doc/ase/parallel.rst`, `doc/tutorials/neb/diffusion.rst`)
- Will all ranks read identical input via `ase.io.read`, or does each rank need a unique file?
- Is backend threading configured (for example Psi4 `num_threads`)? (`doc/ase/calculators/psi4.rst`)
- Is the workload balanced across images/tasks?

### Canonical workflow
1. Confirm MPI-capable environment and communicator discovery (`doc/ase/parallel.rst`).
2. Use world-broadcast-friendly reads when all ranks should share the same atoms (`doc/ase/parallel.rst`).
3. For per-rank files, switch to rank-indexed `Trajectory` reads (`doc/ase/parallel.rst`).
4. For NEB image parallelization, run one image per process as in the tutorial pattern (`doc/tutorials/neb/diffusion.rst`).
5. Validate scaling and load balance with short timing comparisons.

### Minimal working examples
```python
from ase.io import Trajectory
from ase.parallel import world

atoms = Trajectory(f'myfile_{world.rank}.traj')[-1]
```

```bash
mpiexec -np 3 gpaw-python diffusion3.py
```

### Pitfalls
- `ase.io.read` in MPI mode expects all ranks to read the same thing in the same order (`doc/ase/parallel.rst`).
- Per-rank unique-file patterns should use `Trajectory`, not shared `ase.io.read` calls (`doc/ase/parallel.rst`).
- Psi4 defaults to a single thread unless `num_threads` is set (`doc/ase/calculators/psi4.rst`).
- NEB image-parallel runs can still stall if one image is much more expensive than others (`doc/tutorials/neb/diffusion.rst`).
- Restarting without checking trajectory/state consistency can silently waste cycles (`doc/tutorials/neb/diffusion.rst`).

### Convergence/validation checklist
- Verify each rank is assigned intended data/image.
- Check no MPI deadlocks at IO boundaries.
- Compare wall-clock speed versus serial baseline on a small reproducible case.
- Confirm barriers/energies are consistent across restart and parallel reruns.

## Scope
- Handle questions about MPI/OpenMP/GPU execution, scaling, and batch systems.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/calculators/psi4.rst`
- `doc/tutorials/neb/diffusion.rst`
- `doc/ase/parallel.rst`

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
- `ase/calculators/psi4.py`
- `ase/parallel.py`
- `ase/calculators/__init__.py`
- `ase/calculators/vasp/__init__.py`
- `ase/calculators/turbomole/__init__.py`
- `ase/calculators/siesta/__init__.py`
- `ase/calculators/openmx/__init__.py`
- `ase/calculators/lammps/__init__.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
