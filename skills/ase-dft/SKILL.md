---
name: ase-dft
description: This skill should be used when users ask about dft in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Dft

## Scope
- Handle questions about documentation grouped under the 'dft' theme.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/dft/dos.rst`
- `doc/ase/dft/stm.rst`
- `doc/ase/dft/bandgap.rst`

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
- `ase/dft/stm.py`
- `ase/dft/dos.py`
- `ase/dft/bandgap.py`
- `ase/dft/__init__.py`
- `ase/dft/band_structure.py`
- `ase/transport/stm.py`
- `ase/dft/wannierstate.py`
- `ase/dft/wannier.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
