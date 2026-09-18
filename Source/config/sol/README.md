# Local sol2 changes

The bundled headers were generated from sol2 develop commit `c1f95a773c6f8f4fde8ca3efe872e7286afe4444`, which fixes the optional-reference template error exposed by GCC 16.
We also patched associative-container iteration in `sol.hpp` to use the existing `sen()` accessor and initialize the sentinel with `end()` instead of `begin()`.
Preserve these two local fixes when regenerating the headers.
