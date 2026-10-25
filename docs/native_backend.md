# Native perfect-form enumeration

From the checkout root, install the Julia dependencies and start Julia:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=.
```

Compute forms from scratch:

```julia
using ArithmeticVoronoi, Oscar

P = arithmetic_problem(13);
E = enumerate_perfect_forms(P; backend=:native, progress=true);
length(E.representatives)             # 9
E.completion_status                  # :complete
[r.stabilizer_order for r in E.representatives]
E.neighbours                         # labelled facet crossings and transporters
```

`progress=true` prints the rank attained while constructing the starting form,
then the current orbit, facet orbit, number of forms found and number still
awaiting traversal. The total number of orbits is discovered by the computation.
Julia's first import and compilation can take longer than the enumeration.

For presentations, two-cells and matrix-to-word computations, pass the same
backend to the existing pipeline:

```julia
R = run_field_pipeline(12; backend=:native, progress=true,
    directory="results/native_D12");

S = R.simplified[:SL];
F = free_group(length(S.generator_matrices));
G, toG = quo(F, [F(r) for r in S.relators]);
A, toA = maximal_abelian_quotient(G);
elementary_divisors(abelian_invariants(A))  # [2, 6]

P = R.enumeration.problem;
g = matrix(P.field, [1 P.generator; 0 1]);
c = matrix_to_word(R, g; target=:SL);
x = toG(F(Int.(c["word"])));
verify_word_certificate(R, g, c)      # true
```

Use a different output directory when switching backends. Both enumerate the
same mathematical objects, but orbit order, stabilizer generators, facets and
words can have different labels. Certificates refer to their own enumeration.

## Scope and dependencies

The native enumerator handles binary rational forms (`classical_problem()`)
and rank-two forms over maximal quadratic orders on the free lattice
`O^2`. `arithmetic_problem(D)` accepts positive and negative fundamental quadratic
discriminants; negative discriminants use Hermitian forms.
It is not restricted to a table of previously computed fields.
The group-presentation pipeline computes fundamental units and ordinary class
groups with Hecke (`GRH=false`), without a field whitelist. It normalizes the
fundamental unit to be greater than one at the positive square-root embedding;
positive trace and positive coefficient of `w` select it using exact arithmetic.
See the [16-field survey](native_survey.md) for computed examples and timings.
For imaginary fields, the unit group is finite and SL is obtained directly
as the determinant kernel of index 2, 4 or 6. See the
[imaginary quadratic examples](imaginary_quadratic.md) for computations and comparisons.
This is not yet a replacement for
PerfectForms' arbitrary T-spaces or higher ranks.

The computation uses Hecke's exact lattice routines, Nemo arithmetic and
OSCAR's exact polyhedral geometry. It requires no Dutour exporter executable,
CMake build or backend configuration. The PerfectForms C++ adapter and its
unregistered dependency live in `integrations/perfect_forms/`, loaded only
for explicitly requested comparisons. See the [user guide](../USER_GUIDE.md#optional-c-backend).

All enumeration and field-pipeline APIs, including the
[record-based API](requests.md), default to `backend=:native`.
Native enumeration runs afresh,
and has no process deadline or memory limit. `config`, `cache_dir`, `work_dir`,
`timeout` and `force` are C++-backend options. `write_enumeration(path, E)` can
save a native result for inspection; native checkpoint loading is not
implemented. Pipeline output files are still written normally.
The archival `export_dataset`/`verify_dataset` format still requires a raw C++
checkpoint; native enumerations are not yet supported by that archive format.

## Algorithm and checks

The group-presentation and constructive-membership algorithm follows
Braun, Coulangeon, Nebe and Schönnenbeck,
[*Computing in arithmetic groups with Voronoi's algorithm*](https://arxiv.org/abs/1407.6234).
Perfect-form enumeration supplies the vertices and neighbours used by that
algorithm; ridge cycles and stabilizers supply the presentation relations.

The native enumerator is implemented in `src/NativeBackend.jl`. It follows
the relative
Ryshkov-polyhedron approach described by Achill Schürmann in
[Enumerating perfect forms](https://arxiv.org/abs/0901.1587), with Dutour
Sikirić's `polyhedral_common` serving as the independent reference implementation.

1. Start with the trace Gram matrix of the identity field form. Find exact
   shortest vectors and their evaluation matrix on the form space. Deform in
   its nullspace until those evaluations have full rank, then normalize the
   minimum to one.
2. Compute the complete K-linear stabilizer with Hecke, retaining central
   `-I`. Compute the exact cone facets with OSCAR and their stabilizer orbits.
3. For each facet orbit, move in its outward form-space direction to the first
   new minimum. Rational trial points bracket the event; newly short vectors
   give exact bounds on its position. There are no floating-point thresholds.
4. Use determinant and minimal-vector count to filter previously found forms,
   then Hecke's K-linear isometry test to identify the target orbit. New orbits
   enter the traversal queue.
5. Check each common minimal-vector face, outward ray and transporter, cover
   every facet orbit, and verify all reverse incidences before returning
   `completion_status == :complete`.

Forms and matrices pass directly between Julia functions. No JSON, streams or
namelists participate in native enumeration. The returned
`ArithmeticVoronoiEnumeration` has the same mathematical record format as the
validated C++ importer. Transporters satisfy
`H = transpose(g) * F_target * g`; they send neighbour vectors to target
vectors, including their signs. For imaginary fields, replace `transpose`
by conjugate transpose. The rational Gram is the trace of `v* F v`, with
complex conjugation on the first vector; its minimum is normalized to one.

## Tests and C++ comparison

Run the native tests without installing an exporter:

```sh
julia --startup-file=no --project=. test/native.jl
julia --startup-file=no --project=. test/quadratic_units.jl
julia --startup-file=no --project=. test/imaginary.jl
```

These check an explicit exact flip, invalid directions, seed construction,
closed graphs for the classical case and D = 5, 8, 12, 13, stabilizer words,
and a complete D12 GL/SL/PSL pipeline with matrix-to-word certificates and
partial ridge-checkpoint replay. They
deliberately set an invalid exporter configuration to check independence.
The unit tests check all 16 survey fields, including both unit-norm signs,
class numbers one and two, and exact determinant coordinates for powers of
the fundamental unit.

For a fresh comparison against C++, configure an installed exporter and run:

```sh
export PERFECT_FORMS_BACKEND_CONFIG=/absolute/path/to/install/config.json
julia --startup-file=no --project=. test/native_comparison.jl
```

The comparison runs both enumerators for the classical case and D = 5, 8, 12,
13, 17, 21, 24, 28, 29, 33, 37, 40, 41, 44, 56, 60. To select a subset,
append discriminants, for example `test/native_comparison.jl 13 40`.
It matches forms by exact K-linear isometries, compares stabilizer orders and
minimal vectors, then compares every labelled neighbour after transporting
its facet into the reference coordinates. It does not assume identical
representatives or compare counts alone.
