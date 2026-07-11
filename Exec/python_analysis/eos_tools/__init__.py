"""Offline EOS table toolchain for Cerberus tabulated-EOS.

Package layout (see doc/eos_creation_plan.md):
    constants   physical constants (CGS)
    formats     source readers + the frozen .eostab writer/reader
    grids       hull-aware PCHIP regrid onto log-uniform (rho, T)
    condition   derivative blocks, inverse maps, energy shift
    hugoniot    Rankine-Hugoniot locus solver
    qa          QA plots + acceptance checks
    sources     raw-data manifest (fetch/verify)
    cli         command-line entry point (via the eos_table_prep.py shim)

Entry point: Exec/python_analysis/eos_table_prep.py (path-stable shim kept
for the test-suite `run` scripts).
"""
