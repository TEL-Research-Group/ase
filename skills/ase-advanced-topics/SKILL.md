---
name: ase-advanced-topics
description: This skill should be used when users ask about narrow ASE topics that each map to a single documentation page (Atom, Cell, Symbols, Units, Formula, Data, Collections, ULM I/O, MEP module overview, phase diagrams, transport module, or the general module index), with documentation-first routing before source inspection.
---

# ase: Advanced Topics

## Scope
- Handle focused requests that map to one documentation page and do not justify a standalone large-topic skill.
- Route broader follow-ups to neighboring skills (for example modeling, workflows, analysis, or API) when the request expands.

## Route the request
- `atom` basics: `doc/ase/atom.rst`.
- Cell representation and utilities: `doc/ase/cell.rst`.
- Chemical symbols and conversions: `doc/ase/symbols.rst`.
- Unit conventions and CODATA versions: `doc/ase/units.rst`.
- Formula parsing/formatting: `doc/ase/formula.rst`.
- Built-in datasets and references: `doc/ase/data.rst` and `doc/ase/collections.rst`.
- ULM container format details: `doc/ase/io/ulm.rst`.
- MEP module overview: `doc/ase/mep.rst`.
- Phase/Pourbaix diagrams: `doc/ase/phasediagram/phasediagram.rst`.
- Transport module overview: `doc/ase/transport/transport.rst`.
- Top-level module index: `doc/ase/ase.rst`.

## Workflow
- Start with `references/doc_map.md` and cite exact documentation file paths.
- If the user asks for cross-topic workflows (optimization, NEB, MD, calculators), route to the corresponding core skill.
- If docs leave ambiguity, inspect `references/source_map.md` and only then inspect source.
- Prefer targeted source search first (for example: `rg -n "<symbol_or_keyword>" ase`).

## Tutorials and examples
- `doc/tutorials`

## Test references
- `ase/test`

## Optional deeper inspection
- `ase`

## Source entry points for unresolved issues
- `ase/atom.py`
- `ase/cell.py`
- `ase/symbols.py`
- `ase/units.py`
- `ase/formula.py`
- `ase/data/__init__.py`
- `ase/collections/collection.py`
- `ase/io/ulm.py`
- `ase/mep/neb.py`
- `ase/phasediagram.py`
- `ase/transport/calculators.py`
- `ase/__init__.py`
