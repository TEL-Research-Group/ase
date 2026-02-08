---
name: ase-analysis-and-output
description: This skill should be used when users ask about analysis and output in ase; it prioritizes documentation references and then source inspection only for unresolved details.
---

# ase: Analysis and Output

## High-Signal Playbook

### Route the request
- Use this skill for visualization, trajectory handling, vibrational/raman outputs, and Bader-style charge post-processing.
- Route simulation-control questions to `ase-simulation-workflows`.
- Route data-model setup questions to `ase-inputs-and-modeling`.

### Triage questions
- Is the user analyzing trajectories, static snapshots, vibrational spectra, Raman spectra, or charge partitioning? (`doc/ase/io/trajectory.rst`, `doc/ase/visualize/visualize.rst`, `doc/ase/vibrations/vibrations.rst`, `doc/ase/vibrations/raman.rst`, `doc/ase/dft/bader.rst`)
- Is output meant for quick inspection (`ase gui`/`view`) or publication-quality files?
- Are external post-processing tools available (for example Bader executable workflow)?
- Are legacy trajectory formats involved?

### Canonical workflow
1. Choose output artifact and format first (trajectory, image, spectrum table) (`doc/ase/io/trajectory.rst`, `doc/ase/visualize/visualize.rst`).
2. Generate/collect the required intermediate files (vibration data, Raman finite differences, ACF.dat for Bader) (`doc/ase/vibrations/raman.rst`, `doc/ase/dft/bader.rst`).
3. Produce initial plots/visualizations and check units/labels.
4. Convert legacy trajectory data if needed before downstream analysis (`doc/ase/io/trajectory.rst`).
5. Re-run with adjusted sampling/plot windows only after first validated output.

### Minimal working examples
```python
from ase.io.bader import attach_charges

attach_charges(atoms, 'ACF.dat')
for atom in atoms:
    print(atom.symbol, atom.charge)
```

```python
from ase.io import write
write('image.png', atoms)
```

### Pitfalls
- Bader attachment expects compatible charge files (`ACF.dat`) from external tooling (`doc/ase/dft/bader.rst`).
- Molecular vibration outputs include translational/rotational modes that are not true vibrations (`doc/tutorials/tut03_vibrations/vibrations.rst`).
- Old PickleTrajectory files are deprecated and should be converted before reuse (`doc/ase/io/trajectory.rst`).
- Viewer capabilities differ (`ase.gui`, `ngl`, `x3d`, external viewers), so pick tool by context (`doc/ase/visualize/visualize.rst`).
- Raman workflows depend on finite-difference data quality and approximation choices (`doc/ase/vibrations/raman.rst`).

### Convergence/validation checklist
- Confirm generated files open correctly in the chosen viewer/tool.
- Verify units and reference zeroes in plotted quantities.
- Record handling of imaginary/near-zero vibrational modes.
- Ensure analysis scripts are reproducible from saved trajectory/data files.

## Scope
- Handle questions about output formats, analysis, and post-processing.
- Keep responses abstract and architectural for large codebases; avoid exhaustive per-function documentation unless requested.

## Primary documentation references
- `doc/ase/vibrations/raman.rst`
- `doc/ase/visualize/visualize.rst`
- `doc/ase/vibrations/vibrations.rst`
- `doc/ase/dft/bader.rst`

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
- `ase/visualize/plot.py`
- `ase/vibrations/resonant_raman.py`
- `ase/vibrations/raman.py`
- `ase/visualize/__init__.py`
- `ase/vibrations/__init__.py`
- `ase/dft/__init__.py`
- `ase/dft/band_structure.py`
- `ase/visualize/x3d.py`
- Prefer targeted source search (for example: `rg -n "<symbol_or_keyword>" ase`).
