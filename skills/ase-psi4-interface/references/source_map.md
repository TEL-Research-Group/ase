# ase source map: Psi4 Interface

Generated from source roots:
- `ase`

Use this map only after exhausting the topic docs in `references/doc_map.md`.

## Topic query tokens
- `psi4`
- `PSI_SCRATCH`
- `num_threads`
- `method`
- `basis`
- `read`
- `forces`

## Primary implementation
- `ase/calculators/psi4.py`

## Key behaviors to consult in code
- `Psi4.__init__`: imports `psi4` Python module and initializes settings.
- `Psi4.set_psi4`:
  - Handles scratch config (`cfg` / `PSI_SCRATCH`).
  - Applies `reference`, `memory`, `num_threads` (including `"max"` -> CPU count).
  - Rejects unsupported kwargs (`kpts`, `nbands`, `smearing`, `xc`).
  - Maps `method='LDA'` to `svwn`.
  - Builds Psi4 geometry with:
    - `units angstrom`
    - `symmetry <...>`
    - `<charge> <multiplicity>`
    - `no_reorient`
- `Psi4.calculate`:
  - Auto-switches to unrestricted reference when initial magnetic moments exist.
  - Forces path uses `psi4.driver.gradient(...)`; converts units to eV/Angstrom with sign flip.
  - Energy path uses `psi4.energy(...)`; converts Hartree to eV.
  - Writes ASE metadata block (`!ASE Information`) into `label.dat`.
- `Psi4.read`:
  - Reads `label.dat` only if metadata block is present.
  - Restores atoms/parameters/results from embedded JSON.

## Registration and discovery
- `ase/calculators/names.py` (registered calculator name list includes `psi4`).
- `ase/codes.py` (code registry entry for `psi4`).
- `ase/calculators/calculator.py` (`get_calculator_class` special handling for `Psi4`).

## Suggested source entry points
- `ase/calculators/psi4.py`
- `ase/calculators/names.py`
- `ase/codes.py`
- `ase/calculators/calculator.py`
- `ase/test/calculator/psi4/test_psi4_HF_3_21G.py`

## Quick extraction commands
- Command: `rg -n "class Psi4|def set_psi4|def calculate|def read|InputError|PSI_SCRATCH|set_num_threads|driver.gradient|!ASE Information" ase/calculators/psi4.py`
- Command: `rg -n "psi4" ase/calculators/names.py ase/codes.py ase/calculators/calculator.py`
