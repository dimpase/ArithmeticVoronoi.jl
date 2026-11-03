# ArithmeticVoronoi.jl

Compute perfect forms, finite presentations and verified matrix words for
GL(n,O), SL(n,O) and PSL(n,O), n ≥ 2, over Z and maximal quadratic orders
on the free lattice Oⁿ. The default backend uses Julia, OSCAR and Hecke.

Start with the [user guide](guide.md), then see [higher-rank examples](higher_rank.md),
[HAP computations](hap_and_vmh.md), and the [API reference](api.md).

The group algorithm is from Oliver Braun, Renaud Coulangeon, Gabriele Nebe
and Sebastian Schönnenbeck,
[*Computing in arithmetic groups with Voronoi's algorithm*](https://arxiv.org/abs/1407.6234).
The [presentation justification](theory/presentations.md) explains the
implementation and the [validation guides](higher_rank.md) distinguish
completed comparisons from open checks.

Installation does not build C++ or require Magma. Recorded presentations are
available in the independent ArithmeticVoronoiData package and are loaded
only when explicitly requested; see [code and data separation](code_and_data.md).

```@contents
Pages = ["guide.md", "requests.md", "api.md"]
Depth = 2
```
