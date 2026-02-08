---
name: ase-examples-and-tutorials
description: This skill should be used when users ask about examples and tutorials in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Examples and Tutorials

## High-Signal Playbook

### Route the request
- Use this skill when the user needs runnable tutorial patterns and adaptation guidance.
- Route deep theory/method questions to `ase-theory-and-methods`.
- Route backend configuration failures to `ase-calculators`.

### Triage questions
- Which family matches the user problem: NEB/diffusion, defects, QM/MM, vibrations, or calculator-specific examples? (`doc/tutorials/tutorials.rst`, `doc/tutorials/dissociation.rst`, `doc/tutorials/constraints/diffusion.rst`, `doc/tutorials/tut03_vibrations/vibrations.rst`, `doc/tutorials/qmmm/qmmm.rst`)
- Is the request about reproducing a published tutorial result or adapting it to a new system?
- Are required external calculators available for the selected tutorial?
- Is trajectory/barrier post-analysis needed immediately after run?

### Canonical workflow
1. Pick the closest tutorial reference and run it with minimal edits (`doc/tutorials/tutorials.rst`).
2. Confirm expected output artifacts (trajectory, log, barrier/spectrum summary) are produced (`doc/tutorials/dissociation.rst`, `doc/ase/io/trajectory.rst`).
3. Visualize results early (`ase gui`, NEB tools) before parameter sweeps (`doc/tutorials/neb/diffusion.rst`).
4. Modify one parameter family at a time (images, constraints, convergence criteria) and keep a run log.
5. If ambiguity appears, escalate to primary module docs listed in this skill (`doc/ase/phonons.rst`, `doc/ase/vibrations/modes.rst`, calculator docs).

### Minimal working examples
```bash
ase gui neb.traj@-5:
ase nebplot --share-x --share-y neb.traj
```

```python
from ase.io.trajectory import Trajectory

traj = Trajectory('example.traj')
atoms = traj[-1]
```

### Pitfalls
- `doc/tutorials/tutorials.rst` is a legacy index; some content is being moved/ported.
- For simple one-coordinate barriers, NEB may be overkill; constraints-based workflow can be simpler (`doc/tutorials/constraints/diffusion.rst`).
- Skipping visualization can hide path/symmetry issues in NEB tutorials (`doc/tutorials/neb/diffusion.rst`).
- Vibrational exercises include non-physical translational/rotational modes for molecules that must be interpreted correctly (`doc/tutorials/tut03_vibrations/vibrations.rst`).
- Tutorial scripts may depend on external engines and pseudopotential/tool availability.

### Convergence/validation checklist
- Reproduce a baseline tutorial output before adapting.
- Keep trajectory and log files for each variant.
- Validate result trends against tutorial expectations (barrier shape, mode counts, etc.).
- Document all non-default settings used in adapted runs.

## Scope
- Handle questions about worked examples, tutorials, and cookbook usage.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/io/trajectory.rst`
- `doc/ase/calculators/gaussian.rst`
- `doc/tutorials/tutorials.rst`
- `doc/tutorials/defects/defects.rst`
- `doc/ase/calculators/dmol.rst`
- `doc/ase/phonons.rst`
- `doc/ase/vibrations/modes.rst`
- `doc/ase/calculators/nwchem.rst`
- `doc/ase/calculators/gamess_us.rst`
- `doc/ase/calculators/qchem.rst`
- `doc/ase/calculators/gulp.rst`
- `doc/tutorials/dissociation.rst`

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
- `ase/calculators/qchem.py`
- `ase/calculators/nwchem.py`
- `ase/calculators/gulp.py`
- `ase/calculators/gaussian.py`
- `ase/calculators/gamess_us.py`
- `ase/calculators/eam.py`
- `ase/calculators/dmol.py`
- `ase/calculators/castep.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
