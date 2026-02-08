---
name: ase-gui
description: This skill should be used when users ask about gui in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Gui

## Scope
- Handle questions about documentation grouped under the 'gui' theme.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/gui/tools.rst`
- `doc/ase/gui/view.rst`
- `doc/ase/gui/calculate.rst`
- `doc/ase/gui/edit.rst`
- `doc/ase/gui/setup.rst`
- `doc/ase/gui/gui.rst`

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
- `ase/gui/view.py`
- `ase/gui/celleditor.py`
- `ase/gui/atomseditor.py`
- `ase/gui/__init__.py`
- `ase/transport/tools.py`
- `ase/gui/widgets.py`
- `ase/gui/utils.py`
- `ase/gui/ui.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
