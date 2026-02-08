---
name: ase-build-and-install
description: This skill should be used when users ask about build and install in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Build and Install

## Scope
- Handle questions about build, installation, compilation, and environment setup.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/neb.rst`
- `doc/ase/neighborlist.rst`

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
- `ase/gui/po/Makefile`
- `ase/neighborlist.py`
- `ase/dependencies.py`
- `ase/mep/neb.py`
- `ase/cli/build.py`
- `ase/calculators/kim/neighborlist.py`
- `ase/mep/dyneb.py`
- `ase/mep/autoneb.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
