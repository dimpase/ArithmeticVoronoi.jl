# Computation, comparisons and recorded results

The computation engine and recorded-results library are separate Julia
packages in separate Git repositories. Neither depends on the other. The
[data repository](https://github.com/dimpase/ArithmeticVoronoiData.jl) owns the
LFS archives; a current engine checkout has no LFS-managed files. Their boundaries are:

| Location | Purpose | Loaded by the engine? |
| --- | --- | --- |
| `src/` | Arithmetic algorithms, certificate checks and explicit archive I/O | Yes |
| `test/` | Expected values, regression fixtures and corruption tests | No |
| `validation/` and comparison scripts | Published reference presentations and independent comparisons | No |
| [ArithmeticVoronoiData.jl](https://github.com/dimpase/ArithmeticVoronoiData.jl) | Explicitly loaded, versioned precomputed presentations and historical LFS archives | No |
| `results/` | Run outputs, comparison receipts and diagnostic logs | Only user-selected checkpoints or archives, through explicit input |

The engine derives presentation relations, stabilizers, units and homology
inputs from the supplied arithmetic problem. It does not select a known
answer based on the field or rank. Numerical regression expectations and
published presentation tables stay outside `src/`.

## Explicit reuse

Recorded results may be useful inputs to new computations. For example, an
algorithm could use an SL₂ presentation when constructing one for SL₃. The
caller must select and load that record explicitly, then pass it as an input.
Library help must remain optional: every computation accepting it must also
have a complete from-scratch path. That path must work without the library
installed and for inputs absent from its catalogue. Accelerated and unassisted
paths must compute the same mathematical object; their tests should compare
appropriate exact maps or invariants, allowing different generator choices.
An engine entry point accepting such an input should check the properties
its algorithm relies on (group, ring, lattice and any required generator
maps or certificates) and record the library version and source checksum
in its output provenance. It must not silently substitute a library answer
for the requested computation or silently fall back to computation on a
missing library entry.

The current Voronoi GLₙ runner does not use a lower-rank presentation as an
auxiliary input. The data package supplies explicit records for clients and
future algorithms; it does not claim to implement that rank-reduction step.

The current supported reuse mechanisms are:

- `ArithmeticVoronoiData.load_presentation(collection, D; group=...)` loads a
  recorded result, announces the load, and returns source/version/checksum
  provenance. It never performs an arithmetic computation.
- `run_field_pipeline(...; resume=true)` and `run_gl_pipeline(...; resume=true)`
  explicitly replay attaching-loop certificates from the selected run directory.
  Exact geometric validation is required, and replay counts are recorded.
  Presentations and Tietze simplifications are rebuilt from the resulting geometry.
- For the optional C++ enumerator, `force=false` explicitly permits loading
  a checkpoint. The problem, exporter build and complete arithmetic output
  are checked again. Native enumeration always runs from scratch.
- Archive and word-certificate loading requires a supplied path. It is not
  connected to the recorded-results catalogue.

The old automatic archived-baseline selection has been removed. Computation
pipelines do not load saved simplified presentations. The explicit archive
reader checks their recorded matrix maps; that alone is not a new abstract
Tietze proof. Fresh simplification obtains its isomorphism from OSCAR/GAP.

## Regression safeguards

`test/engine_isolation.jl` starts a fresh Julia process with only a copy of
`src/`. There are no reference tools, data package, fixtures or recorded
results in that copy. It computes GL₃(Z) and the real-quadratic D=5 groups,
and checks that deliberately corrupt archived answers do not influence
fresh computation. Explicit replay must reject corrupt loop certificates.

`ArithmeticVoronoiData` has its own tests and CI, without the engine or
OSCAR installed as dependencies. The default engine suite (`Pkg.test()` or
`test/runtests.jl`) uses native computation and a small ordinary-Git fixture;
it needs no data checkout, HAP or C++ executable. HAP runs in a separate CI
job. C++/archive comparisons use the manually dispatched
[integration workflow](../.github/workflows/integration.yml), which pins the
data repository revision. Local bundle checks require explicit paths. Changes to library data
should be reviewed together with independent validation receipts; a test
expectation must never become a special-case answer in the engine.

The [validation receipt](../results/code_data_separation_validation.json)
records 2,623 passing assertions across native, HAP, isolation, data-package
and checkpoint opt-in tests. The real legacy C++ integration check could not
start on this host because its installed binary requires a missing Boost 1.90
shared library. Its separate fixture-based checkpoint test passed; CI builds
the real exporter before running its integration tests.
