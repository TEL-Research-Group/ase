# ASE Psi4 Test Map

## Primary tests
- `ase/test/calculator/psi4/test_psi4_HF_3_21G.py`
  - End-to-end Psi4 calculator usage on `H2O`.
  - Validates:
    - energy magnitude regression,
    - `read('psi4-calc')` restore path,
    - force consistency,
    - analytical vs numerical force agreement.

## Related infrastructure tests and factories
- `ase/test/factories.py`
  - `Psi4Factory` import gating and calculator construction.
- `ase/test/calculator/test_command.py`
  - Shows Psi4 is intentionally excluded from standard command keyword/envvar command tests (uses Python API path and external package dependency).

## Quick extraction commands
- Command: `rg -n -i "psi4|read\\(|forces|numerical" ase/test/calculator/psi4/test_psi4_HF_3_21G.py`
- Command: `rg -n -i "factory\\('psi4'\\)|class Psi4Factory" ase/test/factories.py`
- Command: `rg -n -i "psi4" ase/test/calculator/test_command.py`
