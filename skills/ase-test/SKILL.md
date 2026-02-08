---
name: ase-test
description: This skill should be used when users ask about test in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Test

## Scope
- Handle questions about documentation grouped under the 'test' theme.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `ase/test/testdata/README.md`
- `ase/test/calculator/openmx/md/md_results.txt`

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
- `ase/calculators/openmx/__init__.py`
- `ase/calculators/openmx/writer.py`
- `ase/calculators/openmx/reader.py`
- `ase/calculators/openmx/parameters.py`
- `ase/calculators/openmx/openmx.py`
- `ase/calculators/openmx/dos.py`
- `ase/calculators/openmx/default_settings.py`
- `ase/calculators/calculator.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
