---
name: ase-api-and-scripting
description: This skill should be used when users ask about api and scripting in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: API and Scripting

## High-Signal Playbook

### Route the request
- Use this skill for programmatic API usage: constraints, DB workflows, calculator interfaces, and reusable scripting patterns.
- Route pure execution protocol questions (MD/NEB loops) to `ase-simulation-workflows`.
- Route broad model-construction questions to `ase-inputs-and-modeling`.

### Triage questions
- Is the main task constraint control, database management, or calculator API access? (`doc/ase/constraints.rst`, `doc/ase/db/db.rst`, `doc/ase/calculators/*.rst`)
- Should data be queried from CLI (`ase db`) or Python API? (`doc/ase/db/db.rst`)
- Are filters needed for coupled cell/position optimization? (`doc/ase/filters.rst`)
- Does the workflow need DFT utility objects (band paths/DFT helpers)? (`doc/ase/dft/dft.rst`)

### Canonical workflow
1. Build/load an `Atoms` object and apply constraints explicitly if required (`doc/ase/constraints.rst`).
2. Attach the calculator and run the minimal property call needed (`doc/ase/calculators/emt.rst`, `doc/ase/calculators/orca.rst`, `doc/ase/calculators/turbomole.rst`).
3. Persist results/metadata into ASE DB with clear key-value tags (`doc/ase/db/db.rst`).
4. Query rows with deterministic filters and feed results into analysis/post-processing (`doc/ase/db/db.rst`).
5. Use filter objects when optimizing cell and positions together (`doc/ase/filters.rst`).

### Minimal working examples
```python
from ase.constraints import FixAtoms

c = FixAtoms(mask=[atom.symbol == 'Cu' for atom in atoms])
atoms.set_constraint(c)
```

```python
from ase.db import connect

db = connect('abc.db')
db.write(atoms, relaxed=False)
for row in db.select(relaxed=False):
    print(row.id, row.formula)
```

### Pitfalls
- Constrained atoms remain frozen; direct position edits require clearing constraints first (`doc/ase/constraints.rst`).
- Server-backed DB engines (PostgreSQL/MySQL/MariaDB) require extra backend installation (`doc/ase/db/db.rst`).
- Query expressions differ between exact formula and compositional filters (`doc/ase/db/db.rst`).
- Skipping explicit metadata keys makes later selection brittle (`doc/ase/db/db.rst`).
- Calculator APIs differ by backend; keep interface assumptions local to selected engine docs (`doc/ase/calculators/turbomole.rst`, `doc/ase/calculators/orca.rst`).

### Convergence/validation checklist
- Confirm constraints are applied to the intended atom set.
- Confirm DB rows contain required keys for downstream selection.
- Re-run selection queries and verify deterministic counts.
- Validate that key calculator calls (`energy`, `forces`) are reproducible across script reruns.

## Scope
- Handle questions about language bindings, APIs, and programmatic interfaces.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/constraints.rst`
- `doc/ase/db/db.rst`
- `doc/ase/calculators/turbomole.rst`
- `doc/ase/filters.rst`
- `doc/ase/calculators/emt.rst`
- `doc/ase/utils.rst`
- `doc/ase/calculators/orca.rst`
- `doc/ase/dft/dft.rst`
- `doc/ase/calculators/jacapo.rst`

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
- `ase/calculators/turbomole/__init__.py`
- `ase/calculators/orca.py`
- `ase/calculators/emt.py`
- `ase/calculators/turbomole/writer.py`
- `ase/calculators/turbomole/turbomole.py`
- `ase/calculators/turbomole/reader.py`
- `ase/calculators/turbomole/parameters.py`
- `ase/calculators/turbomole/executor.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
