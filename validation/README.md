# Independent comparisons

This directory contains reference constructions and comparison tools. It is
not loaded by `using ArithmeticVoronoi`, and the engine never calls into it.
Numerical expectations belong in `test/`; recorded presentation data belongs
in the separate [ArithmeticVoronoiData.jl repository](https://github.com/dimpase/ArithmeticVoronoiData.jl). Comparison receipts and
logs belong in `results/`.

To compare a newly computed result with Kirchheimer–Wolfart, explicitly load
the optional comparison module from the repository root:

```julia
include("validation/KirchheimerWolfart.jl")
using .KirchheimerWolfart
compare_kw(R)
```

The published presentations and field-specific reference reductions are
defined only in that module. Generic certificate verification stays in the
engine: it derives its checks from the supplied inputs, without expected
counts, presentations or homology tables.

See [code and data boundaries](../docs/code_and_data.md), including explicit
use of recorded results as auxiliary inputs to future algorithms.
