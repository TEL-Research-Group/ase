---
name: ase-index
description: This skill should be used when users ask how to use ase and the correct topic skill must be selected before going deeper into source code.
---

# ase Skills Index

## Route the request
- Classify the request into one of the available topic skills listed below.
- Prefer abstract, workflow-level guidance for large scientific packages; do not attempt full function-by-function coverage unless explicitly requested.

## Available topic skills
- `ase-getting-started`: Getting Started (initial setup, quickstarts, and core concepts)
- `ase-examples-and-tutorials`: Examples and Tutorials (worked examples, tutorials, and cookbook usage)
- `ase-inputs-and-modeling`: Inputs and Modeling (inputs, system setup, models, and physical parameterization)
- `ase-api-and-scripting`: API and Scripting (language bindings, APIs, and programmatic interfaces)
- `ase-simulation-workflows`: Simulation Workflows (simulation setup, execution flow, and runtime controls)
- `ase-analysis-and-output`: Analysis and Output (output formats, analysis, and post-processing)
- `ase-parallel-hpc`: Parallel and HPC (MPI/OpenMP/GPU execution, scaling, and batch systems)
- `ase-build-and-install`: Build and Install (build, installation, compilation, and environment setup)
- `ase-theory-and-methods`: Theory and Methods (theoretical background and algorithmic methods)
- `ase-calculators`: Calculators (documentation grouped under the 'calculators' theme)
- `ase-psi4-interface`: Focused guidance for ASE+Psi4 setup, limits, API behavior, and troubleshooting
- `ase-gui`: GUI (documentation grouped under the 'gui' theme)
- `ase-vibrations`: Vibrations (documentation grouped under the 'vibrations' theme)
- `ase-spectrum`: Spectrum (documentation grouped under the 'spectrum' theme)
- `ase-dft`: DFT (documentation grouped under the 'dft' theme)
- `ase-test`: Test (documentation grouped under the 'test' theme)
- `ase-advanced-topics`: Consolidated routing for one-doc topics (atom, cell, symbols, units, formula, data, collections, ULM I/O, MEP overview, phase diagrams, transport, and general module index)

## Documentation-first inputs
- `doc/ase`

## Tutorials and examples roots
- `doc/tutorials`

## Test roots for behavior checks
- `ase/test`

## Escalate only when needed
- Start from topic skill primary references.
- If those references are insufficient, search the topic skill `references/doc_map.md`.
- If documentation still leaves ambiguity, open `references/source_map.md` inside the same topic skill and inspect the suggested source entry points.
- Use targeted symbol search while inspecting source (e.g., `rg -n "<symbol_or_keyword>" ase`).

## Source directories for deeper inspection
- `ase`
