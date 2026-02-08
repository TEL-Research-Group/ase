---
name: ase-theory-and-methods
description: This skill should be used when users ask about theory and methods in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Theory and Methods

## Scope
- Handle questions about theoretical background and algorithmic methods.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/eos.rst`
- `doc/ase/dimer.rst`

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
- `ase/eos.py`
- `ase/mep/dimer.py`
- `ase/__init__.py`
- `ase/__main__.py`
- `ase/units.py`
- `ase/thermochemistry.py`
- `ase/symbols.py`
- `ase/stress.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
