# Imaginary quadratic examples

The native backend computes binary Hermitian perfect forms and full
`GL(2,O)`, `SL(2,O)` and `PSL(2,O)` presentations for imaginary quadratic
maximal orders. Its input is the **field discriminant**: use `-4` for
`Q(i)`, `-3` for `Q(sqrt(-3))`, and `-7` for `Q(sqrt(-7))`.

## Compute a group in Julia

From the checkout root, beside `Project.toml`:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=.
```

At the Julia prompt:

```julia
using ArithmeticVoronoi, Oscar

R = run_field_pipeline(-7; backend=:native, progress=true,
    directory="results/my_imaginary_D7");

# Make the Bianchi group PSL(2,O) as an OSCAR finitely presented group.
S = R.simplified[:PSL];
F = free_group(length(S.generator_matrices));
G, toG = quo(F, [F(r) for r in S.relators]);
gens(G)
relators(G)
S.generator_matrices                 # exact determinant-one lifts

A, toA = maximal_abelian_quotient(G);
elementary_divisors(abelian_invariants(A))   # [2, 0]: C2 × Z

# Express a matrix as a word in precisely these generators.
P = R.enumeration.problem;
g = matrix(P.field, [1 P.generator; 0 1]);
c = matrix_to_word(R, g; target=:PSL);
x = toG(F(Int.(c["word"])));
verify_word_certificate(R, g, c)      # true
```

Select `R.simplified[:SL]` or `R.simplified[:GL]` and use the matching
`target=:SL` or `target=:GL` to work with those groups. GL and SL keep `-I`.
PSL is `SL/{I,-I}`, so its matrix words are understood modulo sign.

For perfect forms alone:

```julia
E = enumerate_perfect_forms(arithmetic_problem(-23);
    backend=:native, progress=true);
length(E.representatives)            # 3
[r.stabilizer_order for r in E.representatives]
```

No C++ build or backend configuration is needed. The first Julia import and
compilation can take minutes; subsequent computations are much quicker.
The runner computes fresh results by default. Only `resume=true` enables
validated replay of saved ridge certificates from the selected directory.
Native enumeration and presentation simplification always run afresh.
Use a new directory for a fully fresh computation.

## Five-field survey

All five fields were computed from scratch, including ridge cycles, full
presentations, Tietze simplifications and exact matrix-to-word certificates.
Here `Cn` means a cyclic group of order `n`, and the last three columns are
abelianizations, not the original infinite groups.

| Discriminant | Class number | Perfect forms | GL abelianization | SL abelianization | PSL abelianization |
| ---: | ---: | ---: | :--- | :--- | :--- |
| -3 | 1 | 1 | C6 | C3 | C3 |
| -4 | 1 | 1 | C2 × C4 | C2 × C2 | C2 × C2 |
| -7 | 1 | 1 | C2³ | Z × C4 | Z × C2 |
| -15 | 2 | 2 | C2⁴ | Z² × C12 | Z² × C6 |
| -23 | 3 | 3 | C2⁵ | Z³ × C12 | Z³ × C6 |

To reproduce the survey and its comparisons:

```sh
julia --threads=1 --startup-file=no --project=. scripts/run_imaginary_fields.jl
```

This writes full certificates under `results/imaginary_fields/D*/`, a
`summary.json`, per-field `comparison.json` files, and a compact
`presentations.json`. It compares the polytope types with Yasaki and SL/PSL
abelianizations with HAP's independent complexes through Oscar/GAP.
HAP is needed only by this comparison script, not by the native group pipeline.
Use `--output=DIR` to change the destination, `--resume` to retain completed
fields, or append a subset such as `-4 -23`.

The committed [presentation bundle](https://github.com/dimpase/ArithmeticVoronoiData.jl/blob/main/data/imaginary_presentations.json)
contains all 15 simplified presentations with exact matrices, relators,
comparison receipts and timings. Full local certificates for this survey
are in `results/imaginary_survey_20261021/`; the bundle does not embed them.

Warm native enumeration times (median of three fresh runs after a warm-up)
were:

| Discriminant | Perfect-form enumeration, seconds |
| ---: | ---: |
| -3 | 0.0290 |
| -4 | 0.0244 |
| -7 | 0.0183 |
| -15 | 0.0370 |
| -23 | 0.0510 |

These timings cover perfect forms, stabilizers and labelled neighbours,
excluding group presentations and Julia startup/compilation. The machine
was an Intel Core Ultra 7 165U with Julia 1.13.2 and OSCAR 1.8.2, using one
Julia thread and one BLAS thread. The bundle also records the latest full
pipeline stage times; those reruns replayed existing ridge and simplification
certificates and are not fresh end-to-end benchmarks.

The batch checked 147 matrix-to-word certificates. The targeted model,
unit, real/imaginary native pipeline and portable-bundle tests passed
1,773 assertions in total.

## Independent checks and Page's work

The perfect-form counts **and polytope types** agree with
[Yasaki, Tables 2–4](https://arxiv.org/abs/0908.1762): a tetrahedron for -3,
an octahedron for -4, a triangular prism for -7, a tetrahedron and an
octahedron for -15, and an octahedron, triangular prism and square pyramid
for -23. Unit multiples of minimal vectors represent the same polytope vertex.

All ten SL/PSL abelianizations agree with the independent complexes in
[HAP](https://gap-packages.github.io/hap/tutorial/chap14.html), using
`ResolutionSL2QuadraticIntegers` and `ResolutionPSL2QuadraticIntegers`.
This checks abelianizations; it does not claim a certified isomorphism
between our presentations and HAP's presentations.

[Aurel Page's paper](https://www.normalesup.org/~page/Recherche/Documents/articles/kln_gps.pdf),
Section 3.2.1 and Table 3.1, includes three of these fields. The historical
KleinianGroups timings for PSL/PGL are respectively 0.93/0.83 seconds for
-3, 0.92/2.32 for -15, and 1.22/2.06 for -23. These are on different hardware
and measure a different computation; they do not establish a speed ratio.
PGL is GL modulo **all** scalar units, and is not PSL.

Page's unmodified KleinianGroups 1.0 was rerun with Magma 2.30-1 on this
machine. All 15 GL/SL/PSL abelianizations for these five fields agree with
the native results above. The separate
[Page presentation bundle](https://github.com/dimpase/ArithmeticVoronoiData.jl/blob/main/data/page_presentations.json) contains
the exact generator matrices, relators, abelian invariants and CPU timings.
The Page entries in the original native survey bundle record the historical
paper timings, not this rerun.

The comparison driver, [scripts/compare_page.m](../scripts/compare_page.m),
uses the entry points in Page's own `example.m`: `BianchiOrder`,
`NormalizedBasis`, `Presentation`, `LiftPresentation` and `QuatToMatrix`.
Each independent computation constructs a fresh order and resets the random
seed to 1 for reproducibility. Clone the data repository beside this engine
checkout as shown in the [test instructions](../USER_GUIDE.md#tests) to supply
the explicitly selected native comparison bundle.
The driver prints the results and exports JSON in the same matrix and word
encoding as our native bundle. With a licensed Magma installation, download the
[unmodified KleinianGroups 1.0 archive](https://www.normalesup.org/~page/Recherche/Logiciels/KleinianGroups/KleinianGroups-1.0.tar.gz)
and run from its parent directory:

```sh
# Run these lines from this repository's root.
AV_REPO="$PWD"
mkdir -p /tmp/av-page
cd /tmp/av-page
curl -fLO https://www.normalesup.org/~page/Recherche/Logiciels/KleinianGroups/KleinianGroups-1.0.tar.gz
tar -xzf KleinianGroups-1.0.tar.gz
magma -b "OutputFile:=/tmp/av-page/fresh_page_presentations.json" "$AV_REPO/scripts/compare_page.m"
cd "$AV_REPO"
julia --threads=1 --startup-file=no --project=. test/page_comparison.jl /tmp/av-page/fresh_page_presentations.json ../ArithmeticVoronoiData.jl/data/imaginary_presentations.json
```

Page's `BianchiOrder(d)` takes the positive radicand in `Q(sqrt(-d))`:
for discriminant -4 the script passes `1`.
The Julia check verifies all exported matrix relations exactly (modulo
sign for PSL), recomputes the abelianizations with OSCAR and compares them
with the native bundle. These checks do not prove an isomorphism between
the two finite presentations.
The recorded run passed all 93 checks, covering 250 generator matrices and
512 relators. Page's domain and presentation computations took 1.1–2.2 CPU
seconds per group on this host, including the lift for SL/GL and excluding
order construction, export and abelianization. These measure more work than
the native enumeration-only timings above.

## Load a saved presentation

To use the compact bundle without recomputing its geometry, install the data
package in a client environment as shown in the [user guide](../USER_GUIDE.md#load-saved-presentations),
then run:

```julia
using ArithmeticVoronoiData, Oscar

entry = load_presentation(:native_imaginary, -23; group=:PSL);
F = free_group(length(entry.generator_matrices));
G, toG = quo(F, [F(Int.(r)) for r in entry.relators]);
A, toA = maximal_abelian_quotient(G);
elementary_divisors(abelian_invariants(A))   # [6, 0, 0, 0]
```

Matrix entries in the bundle are pairs of rational strings `[a,b]`, meaning
`a+b*w` with `w=(s+sqrt(D))/2` and `s=D mod 4`. A relator is a list of signed
one-based generator indices, evaluated in written matrix-product order.
Geometric matrix-to-word reduction needs the full `R` from a pipeline run.

## Implementation and tests

The imaginary form space has rational dimension four. The trace Gram is
`Tr(v† F v)` on the basis `(e1,w*e1,e2,w*e2)`, and the group
acts by `F -> g† F g`, where `†` denotes conjugate transpose. Hecke supplies exact Hermitian lattice
isometries and automorphism groups; the existing native facet traversal and
ridge machinery operate on the rational trace forms.

The unit group has order 6 for -3, 4 for -4 and 2 otherwise. SL extraction
uses the full finite determinant kernel, with cyclic Reidemeister–Schreier
rewriting. There is no quotient by an infinite fundamental scalar in this
case. PSL is formed only after SL and its central `-I` have been verified.

```sh
julia --startup-file=no --project=. test/imaginary.jl
julia --startup-file=no --project=. test/imaginary_bundle.jl ../ArithmeticVoronoiData.jl/data/imaginary_presentations.json
julia --startup-file=no --project=. test/native.jl
julia --startup-file=no --project=. test/quadratic_units.jl
```

The imaginary tests check Hermitian congruence, exclude semilinear
conjugation from the K-linear group, verify the published polytope types,
and run all five GL/SL/PSL pipelines and both raw and simplified word maps.
The bundle test loads all 15 saved presentations into OSCAR, checks exact
matrix relations and recomputes their abelianizations.
