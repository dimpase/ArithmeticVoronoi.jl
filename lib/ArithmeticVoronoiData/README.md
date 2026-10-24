# ArithmeticVoronoiData

An explicit library of **precomputed** arithmetic-group presentations. This
is a separate Julia package, maintained in this repository for now. It has
no dependency on ArithmeticVoronoi, OSCAR, HAP or Magma, and the computation
engine has no dependency on it. It is not registered yet.

From the repository root, install it in the Julia environment where you want
to use recorded results:

```julia
using Pkg
Pkg.activate(; temp=true)  # a client environment, separate from either package
Pkg.develop(path="lib/ArithmeticVoronoiData")
using ArithmeticVoronoiData
available_presentations(:native_imaginary)
R = load_presentation(:native_imaginary, -7; group=:PSL);
R.origin                        # :recorded_library
R.provenance                    # source file, checksum and library version
```

The collections are `:native_real` (48 GL₂/SL₂/PSL₂ presentations over 16
real quadratic fields), `:native_imaginary` (15 over five imaginary quadratic
fields), and `:page` (15 independently computed with Aurel Page's Magma code).
Every load reports that it is reading recorded data. A missing entry is an
error; loading data never falls back to computation.

`R.generator_matrices` contains rational coefficient strings `[a,b]` for
`a+b*w`, where `w=(s+sqrt(D))/2` and `s=D mod 4`. Relators are signed, one-based
generator indices. PSL matrices are SL lifts, with relators evaluating to ±I.
Use OSCAR separately to construct a finitely presented group:

```julia
Pkg.add(PackageSpec(name="Oscar", version="1.8.2"))
using Oscar
F = free_group(length(R.generator_matrices));
G, toG = quo(F, [F(w) for w in R.relators]);
```

To compute a new answer, use ArithmeticVoronoi's `run_field_pipeline` or
`run_gl_pipeline`. Those functions do not consult this library.
Recorded presentations can also be supplied explicitly to algorithms that
accept them as auxiliary inputs. Such algorithms should retain their
provenance and check the required arithmetic compatibility. The current
GLₙ runner does not yet use lower-rank presentations in this way.
Such use must remain optional: a computation must also be possible from
scratch, with no library installed and no matching recorded entry.

The files under `data/` retain their original schemas and provenance. Tests
and validation reports belong to the computation repository, not this catalogue.
