---
name: ase-vibrations
description: This skill should be used when users ask about vibrations in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Vibrations

## Scope
- Handle questions about documentation grouped under the 'vibrations' theme.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/vibrations/franck_condon.rst`
- `doc/ase/vibrations/H2O_EMT_summary.txt`
- `doc/ase/vibrations/infrared.rst`

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
- `ase/vibrations/franck_condon.py`
- `ase/vibrations/infrared.py`
- `ase/vibrations/__init__.py`
- `ase/vibrations/vibrations.py`
- `ase/vibrations/resonant_raman.py`
- `ase/vibrations/raman.py`
- `ase/vibrations/placzek.py`
- `ase/vibrations/pickle2json.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
