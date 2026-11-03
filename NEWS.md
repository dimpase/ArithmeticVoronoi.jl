# Release notes

## 0.1.0

First release proposed for Julia's General registry. Computes perfect forms,
GL(n,O), SL(n,O) and PSL(n,O) presentations and matrix words for n ≥ 2,
on free lattices over Z and maximal quadratic orders. PSL uses the full
scalar centre, including nonreal roots of unity in higher rank.
Includes a native provider for HAP cell complexes.

Installation uses registered Julia dependencies. Computation defaults to the
native Julia/Hecke/OSCAR backend, with output under the caller's working
directory. The C++ comparison adapter uses a separate, explicitly loaded
environment. Recorded presentations remain in the optional
ArithmeticVoronoiData package and are never loaded implicitly by the engine.

Registration review updates: compatibility bounds allow Julia 1.10 and later,
OSCAR 1.8 and later within 1.x, Hecke 0.39, and JSON 1.x. JSON replaces the
deprecated JSON3 dependency. The root manifest is no longer distributed,
and subprocess tests use the active client test environment. CI tests LTS
and stable Julia with coverage. Certificate hashing uses an owned byte buffer
to avoid repeated string hashing on Julia 1.10, with unchanged SHA-256 digests. The manual, algorithm and audit notes now
build with Documenter, including exported API docstrings.
