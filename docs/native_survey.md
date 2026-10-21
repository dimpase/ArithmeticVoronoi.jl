# Real quadratic survey

The survey uses field discriminants
`5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 37, 40, 41, 44, 56, 60`.
The order is the maximal order `O` and the lattice is the free module `O^2`.
For example, discriminant 60 means `Q(sqrt(15))`.

## Run the full computation

From the checkout root:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --threads=1 --startup-file=no --project=. scripts/run_native_fields.jl
```

This computes all 16 fields using Julia, Hecke and OSCAR. Each field gets
perfect forms, complete ridge coverage, GL/SL/PSL presentations, certified
Tietze simplifications, abelianizations and exact matrix-to-word certificates.
The word checks include upper and lower unipotents, unit diagonals, scalars,
`-I` and a product of these matrices, in raw and simplified alphabets.
No exporter configuration or C++ build is needed.

The script prints stages and ridge progress and updates
`results/native_fields/summary.json` after each field. Its `status` is
`complete` only when all requested fields and word checks pass. An individual
failure is recorded while the remaining fields continue. The process exits
unsuccessfully if any field failed.

Select a subset or a different destination:

```sh
julia --threads=1 --startup-file=no --project=. scripts/run_native_fields.jl \
    --output=results/my_native_batch 17 21 60
```

Each `D*/` directory contains `field.json`, `enumeration.json`,
`GL_raw.json`, `SL_raw.json`, `PSL_raw.json`, the three `*_simplified.json`
files, `two_cells/`, `words/`, `word_validation.json` and `run.json`.
Simplified presentation files contain the exact generator matrices, relators
and both directions of the generator maps. Results directories are ignored
by Git; the script regenerates them. A rerun recomputes native perfect forms
and can reuse compatible ridge and simplification certificates in its directory.
Add `--resume` to retain fields already marked complete in that batch's
`summary.json` and continue the unfinished fields. Loop certificates are saved
individually; the aggregate coverage file is refreshed about once per second
and always on completion, a requested stop or an error.

## Use the computed groups in Julia

For a live result with matrix-to-word operations:

```julia
using ArithmeticVoronoi, Oscar

R = run_field_pipeline(60; backend=:native,
    directory="results/native_fields/D60", reuse_baseline=false);
S = R.simplified[:SL];
F = free_group(length(S.generator_matrices));
G, toG = quo(F, [F(r) for r in S.relators]);

P = R.enumeration.problem;
g = matrix(P.field, [1 P.generator; 0 1]);
c = matrix_to_word(R, g; target=:SL);
x = toG(F(Int.(c["word"])));
verify_word_certificate(R, g, c)  # true
```

To load just an already computed finitely presented group:

```julia
using Oscar, JSON3

data = JSON3.read(read("results/native_fields/D60/SL_simplified.json", String));
F = free_group(length(data.generators));
G, toG = quo(F, [F(Int.(r)) for r in data.relators]);
A, toA = maximal_abelian_quotient(G);
elementary_divisors(abelian_invariants(A))
```

Use `GL` or `PSL` in place of `SL` to select the other presentations. PSL
generator matrices are lifts to SL, and its matrix relations hold modulo
`{I,-I}`. Generator indices always refer to the same presentation file.

## Enumeration results, 21 October 2026

All 16 native enumerations agree with fresh C++ computations: 898 perfect-form
orbits and 5,667 labelled neighbours in total. The comparison matches every
form by exact K-linear isometry, checks stabilizers and minimal vectors, and
matches every facet crossing after changing coordinates. It does not rely on
counts alone.

| Discriminant | Orbits | Labelled neighbours | Native (s) | C++ process (s) | C++ full API (s) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 5 | 2 | 8 | 0.035 | 0.063 | 0.124 |
| 8 | 2 | 7 | 0.040 | 0.063 | 0.325 |
| 12 | 3 | 15 | 0.064 | 0.105 | 0.939 |
| 13 | 9 | 46 | 0.143 | 0.273 | 0.417 |
| 17 | 34 | 200 | 0.516 | 0.988 | 1.519 |
| 21 | 16 | 90 | 0.248 | 0.527 | 0.825 |
| 24 | 22 | 135 | 0.347 | 0.673 | 1.005 |
| 28 | 30 | 200 | 0.547 | 1.116 | 1.754 |
| 29 | 39 | 226 | 0.586 | 1.262 | 1.847 |
| 33 | 149 | 876 | 2.394 | 4.818 | 9.734 |
| 37 | 97 | 554 | 1.417 | 3.032 | 4.979 |
| 40 | 52 | 372 | 1.077 | 2.105 | 3.056 |
| 41 | 144 | 932 | 2.557 | 5.877 | 11.052 |
| 44 | 66 | 465 | 1.379 | 2.674 | 3.947 |
| 56 | 144 | 914 | 2.257 | 4.924 | 9.882 |
| 60 | 89 | 627 | 1.626 | 3.600 | 5.414 |

Native times are medians of three fresh enumerations after a warm-up for each
field, with garbage collection before each measurement and progress disabled.
C++ times are one fresh run per field; the full API additionally imports,
validates and saves the export. Neither backend reuses saved enumeration
results. The C++ bridge was warmed on D5; these single-run full-API measurements
can still include compilation of validation paths specific to later fields.

Environment: Intel Core Ultra 7 165U, Julia 1.13.2, OSCAR 1.8.2, one Julia
thread and one BLAS thread. The existing C++ exporter uses GNU 16.2.0 with
`-O1 -DNDEBUG`; these are not comparisons against an `-O3` rebuild. Importing
the Julia packages took 13.59 seconds; the first native call took 40.37 seconds,
including compilation. Those costs are excluded from the warm table.

The local receipts, exact native enumerations, C++ checkpoints and benchmark
driver are in `results/discriminant_survey_20261021/`. To repeat the exact
comparison with an installed exporter:

```sh
export PERFECT_FORMS_BACKEND_CONFIG=/absolute/path/to/install/config.json
julia --startup-file=no --project=. test/native_comparison.jl
```

The comparison also includes the classical binary rational case. Append
discriminants to run a subset. Fundamental units, determinant coordinates and
ordinary class groups are tested separately by `test/quadratic_units.jl`.
The ordinary class number is two for D40 and D60 and one for the other fields
in this survey; Hecke computes these with `GRH=false`.

## Full group presentations

All 48 GL/SL/PSL presentations were constructed and simplified with certified
Tietze maps. The batch also verified 432 matrix-to-word certificates, and
all 48 saved simplified presentations were loaded directly into OSCAR and
their abelianizations recomputed.

Each size below is **generators / relators**. These are the sizes obtained
in this run; different valid GAP simplifications may have different sizes.

| Discriminant | GL | SL | PSL |
| ---: | ---: | ---: | ---: |
| 5 | 2 / 12 | 2 / 6 | 2 / 6 |
| 8 | 3 / 10 | 2 / 5 | 2 / 5 |
| 12 | 3 / 13 | 2 / 6 | 2 / 6 |
| 13 | 3 / 12 | 2 / 5 | 2 / 7 |
| 17 | 4 / 15 | 2 / 8 | 2 / 6 |
| 21 | 3 / 10 | 3 / 9 | 3 / 8 |
| 24 | 4 / 21 | 3 / 15 | 3 / 11 |
| 28 | 4 / 18 | 3 / 13 | 2 / 9 |
| 29 | 3 / 18 | 2 / 17 | 2 / 15 |
| 33 | 4 / 19 | 3 / 18 | 3 / 18 |
| 37 | 3 / 18 | 3 / 15 | 3 / 13 |
| 40 | 3 / 17 | 2 / 17 | 3 / 15 |
| 41 | 4 / 24 | 2 / 19 | 3 / 12 |
| 44 | 3 / 20 | 2 / 21 | 2 / 15 |
| 56 | 4 / 27 | 2 / 24 | 3 / 20 |
| 60 | 3 / 34 | 2 / 23 | 2 / 20 |

`Cn` denotes a cyclic group of order `n`; `1` denotes the trivial group.

| Discriminant | GL abelianization | SL abelianization | PSL abelianization |
| ---: | :--- | :--- | :--- |
| 5 | Z × C2 | 1 | 1 |
| 8 | Z × C2 × C2 | C2 × C2 | C2 × C2 |
| 12 | Z × C2 × C2 | C2 × C6 | C2 × C6 |
| 13 | Z × C2 | C3 × C3 | C3 × C3 |
| 17 | Z × C2 × C2 × C2 | C4 × C4 | C2 × C4 |
| 21 | Z × C2 | C3 | C3 |
| 24 | Z × C2 × C2 × C2 | C2 × C6 | C2 × C6 |
| 28 | Z × C2 × C2 | C6 × C6 | C6 × C6 |
| 29 | Z × C2 | 1 | 1 |
| 33 | Z × C2 × C2 × C2 | C4 × C12 | C2 × C12 |
| 37 | Z × C2 | C3 × C3 | C3 × C3 |
| 40 | Z × C2 × C2 | C6 × C6 | C6 × C6 |
| 41 | Z × C2 × C2 × C2 | C4 × C4 | C2 × C4 |
| 44 | Z × C2 × C2 | C2 × C2 | C2 × C2 |
| 56 | Z × C2 × C2 × C2 | C2 × C2 | C2 × C2 |
| 60 | Z × C2 × C2 | C2 × C6 | C2 × C6 |

The compact [presentation bundle](../results/native_presentations.json)
contains every simplified generator matrix and relator, the exact field
definitions, abelianizations and enumeration timings. To load one of those
groups without rerunning the geometric computation:

```julia
using Oscar, JSON3

bundle = JSON3.read(read("results/native_presentations.json", String));
data = bundle.fields["60"].groups.SL;
F = free_group(length(data.generators));
G, toG = quo(F, [F(Int.(r)) for r in data.relators]);
```

The full local run directories are `results/native_fields_20261021/D*/`.
They also retain raw presentations, both directions of the Tietze maps,
geometric certificates and the word certificates. The compact bundle does
not include all of that proof data and does not use the older C++ archive
format accepted by `verify_dataset`.

Per-stage group timings are in the batch `summary.json`. D5 includes first-use
compilation. The first six fields used the original per-ridge ledger writing;
the remaining ten used periodic ledger updates. D24 also includes compilation
after restarting Julia and reused its completed loop files. These are run
receipts, not uniform fresh-run group benchmarks.

A separate fresh D17 check reduced the two-cell stage from 178.52
to 23.15 seconds after the checkpoint change, while producing
all 1,352 identical loop certificates. The native test suite also checks that
explicit stops persist coverage and that a resumed run replays saved loops.
