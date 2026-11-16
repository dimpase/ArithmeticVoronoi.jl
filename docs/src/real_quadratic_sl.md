# Certified real-quadratic SLₙ complexes

The determinant-height construction accepts every rank `n ≥ 2` over a real
quadratic maximal order O, on the free lattice Oⁿ. It constructs an
SLₙ(O)-equivariant regular level of the GLₙ Voronoi complex and certifies its
contractibility by exact ascending and descending link collapses. The full
level has dimension `n²−1`. Every run must pass certification; an accepted rank
is not a promise that enumeration or certification will finish quickly.

Rank two was introduced in ArithmeticVoronoi 0.1.1. Higher ranks require the
subsequent development version. The completed end-to-end regression cases are
SL₂ over Q(√5) and Q(√3), with discriminants 5 and 12; see the validation section
for the separate higher-rank checks.

## Run a fresh computation

For the small rank-two example:

```julia
using ArithmeticVoronoi, Oscar
import Oscar: GAP
K = hap_voronoi_complex(5; group=:SL, resolution_length=4, progress=true);
R = GAP.Globals.FreeGResolution(K, 4);
C = GAP.Globals.TensorWithIntegers(R);
[GAP.gap_to_julia(GAP.Globals.Homology(C,n)) for n in 1:3]
```

The registered HAP provider accepts the same request, with a squarefree
radicand in HAP's convenience call:

```julia
enable_hap_julia_backend(; resolution_length=4)
K = GAP.evalstr("VoronoiComplexSL(2,5 : backend := \"julia\")");
```

To inspect the geometric cells, reduced chain model and exact witnesses:

```julia
result = ArithmeticVoronoi.RealQuadraticSL.compute_real_sl_complex(5; rank=2);
length.(result.geometric.cells)
length.(result.chain.cells)
result.certificate
```

The original command-line runner remains available in a source checkout:

```sh
julia --project=. experimental/RealQuadraticSL/run.jl 5 /absolute/path/to/output 2
```

Set `ARITHMETIC_VORONOI_HAP_ROOT` to a HAP checkout to select it explicitly.
For constructing the complex, an installed HAP also suffices. No Magma or C++
is called. Enumeration and certification run afresh. The optional output directory
receives link data, exact collapse witnesses and construction timings; these files
are outputs only and are never loaded implicitly. Use `reduce=false` to export
the full geometric complex to HAP instead of its smaller chain-equivalent model.

For rank three, the corresponding calls are:

```julia
K = hap_voronoi_complex(5; rank=3, group=:SL, resolution_length=2, progress=true);
# Or, after enable_hap_julia_backend(; resolution_length=2):
K = GAP.evalstr("VoronoiComplexSL(3,5 : backend := \"julia\")");
```

`resolution_length=m` computes only the cells through dimension `m` when the
full level is larger. It still certifies **all** ascending and descending links
of the full level. The returned HAP object supports a free resolution through
that degree, hence homology through degree `m−1`. Its
`arithmeticVoronoi.computedCellDimension` and `certifiedCellDimension` metadata
record the distinction. Asking it for uncomputed higher cells raises an error.
A later call with a larger resolution length computes a larger skeleton afresh.

The detailed `compute_real_sl_complex` interface constructs the full level by
default. Set `max_cell_dimension=m` explicitly to request a skeleton; then
`resolution_length` must not exceed `m` if the skeleton is proper. A skeleton
is not itself asserted to be contractible. `reduce=false` keeps its geometric
cells instead of performing algebraic chain cancellations.

The level construction requires a positive quadratic field discriminant and the
free lattice over its maximal order. It returns no complex if its chosen
height cannot be certified. This is not a proof that such a complex does not exist
for that field. A contracting homotopy for HAP cup-product calculations is not
provided.

## Construction and certificate

Let W be the bounded Ryshkov cell complex already computed for GLₙ(O), and let
ε generate the infinite part of O×. The existing well-rounded-retract theorem
makes W contractible, with finitely many cell orbits and finite stabilizers.
Pass to its barycentric subdivision: a simplex is a flag of cells of W.

Write `χ(g)` for the exponent of ε in `det(g)`. For each cell-orbit representative
c choose a rational number b(c), and define the height of the corresponding
barycentric vertex by

```
h(g c) = b(c) + χ(g).
```

This is well defined: a finite cell stabilizer has determinant exponent zero,
which the implementation checks. Different orbit representatives receive distinct
fractional parts, so every edge has distinct endpoint heights. Extend h affinely
on each flag simplex. Its vertex heights form a discrete subset of the rationals,
and GLₙ(O) acts by integer translations on height.

The rational b(c) are proposed from approximate determinant heights of invariant
interior forms, with small rational perturbations. Approximation is used only to
choose a candidate. All subsequent comparisons, unit exponents and certificates
are exact.

The construction enumerates the entire incidence poset at each barycentric
vertex orbit. Comparable cells form edges; chains form simplices, so each link
is the flag complex of its graph. Its ascending and descending subcomplexes
are induced flag complexes. This avoids expanding every maximal flag.

The certificate records dominated-vertex and dominated-edge removals. Each
removal has a cone link, checked by exact neighbour-set inclusion, and expands
into elementary simplicial collapses. An independent replay checks every cone
condition and requires exactly one vertex to remain. See
[Boissonnat–Pritam, *Edge Collapse and Persistence of Flag Complexes*](https://doi.org/10.4230/LIPIcs.SoCG.2020.19).
No numerical topology test is substituted for the certificate. Failure to
produce a witness aborts construction.

By the Morse lemma and Corollary 2.6 of
[Bestvina–Brady, *Morse theory and finiteness properties of groups*, §2](https://people.math.osu.edu/davis.12/courses/8800-20/Bestvina-Brady.pdf),
these link certificates imply that a regular level has the homotopy type of W.
Thus `h⁻¹(0)` is contractible. SLₙ(O) preserves it. The height kernel has
only the additional determinant sign, so SLₙ has finitely many cell orbits
in this level.

The intersection of a mixed-height simplex with the level is a product of two
simplices. Its facets are the mixed-height subflags obtained by deleting one
vertex. Flags have no inversions, hence their stabilizers fix the level cells
pointwise. Exact determinant-one orbit splitting supplies the SL transporters
and stabilizers. The real scalar centre is `{I,−I}` for even `n` and `{I}`
for odd `n`; −I is excluded from SL in odd rank. Boundary signs use the height
coorientation, with zero-cells normalized to positive augmentation.

The implementation verifies equivariant `∂²=0`, augmentation, and invariance of
each boundary under its stabilizer. Optional algebraic cancellations remove pairs
of isomorphic permutation modules; conjugacy of the stabilizers is checked at
every pivot. The full geometric object is retained separately. The reduced
object is a chain model, not a claim that its basis still consists of geometric
product cells. Both can be exported as HAP non-free resolutions.

## Validation and separation from computation

The HAP integration CI tests `hap_voronoi_complex` and the registered GAP
provider on discriminants 5 and 12, checking integral homology through degree
three. These tests are in `test/hap_real_quadratic_sl.jl`.

Extended optional tests are in `test/real_quadratic_sl.jl`; they are not in the
default small test suite. They check certificate replay, reject incomplete
witnesses and changed boundary orientations, and compare H₁ with an
independently constructed SL presentation.

The independent comparison in rank `n` restricts a GLₙ complex to
`det⁻¹(⟨εⁿ⟩) = SLₙ(O) × ⟨εI⟩`, an index-`2n` subgroup. The integral Künneth
formula gives `H_k(SLₙ × Z) = H_k(SLₙ) ⊕ H_(k−1)(SLₙ)`. The subgroup of square
determinants gives this product only for `n=2`.

The same comparison can use freshly generated Magma/VMH GLₙ complexes, with
explicit conversion of VMH's rational row representation. Magma does not
construct the real-quadratic SLₙ level in this comparison. The reference runner
accepts `(D,rank)` cases; computation remains independent of Julia's cells.

```sh
# Optional: licensed Magma plus a VMH checkout; output must be new.
julia --project=. test/prepare_real_sl_magma.jl /path/to/VMH-TotallyRealNumberFields /path/to/references

# Uses HAP with RestrictNonFreeResolution; no precomputed expected answers.
ARITHMETIC_VORONOI_HAP_ROOT=/path/to/hap julia --project=. test/real_quadratic_sl.jl /path/to/test-output /path/to/references
```

For a higher-rank reference, call the runner's Julia function explicitly:

```julia
include("test/prepare_real_sl_magma.jl")
prepare_real_sl_magma("/path/to/VMH-TotallyRealNumberFields", "/path/to/new-references";
                     cases=((D=5,rank=3),), deduplicate_faces=true)
include("test/real_quadratic_sl.jl")
run_real_sl_tests("/path/to/new-test-output"; magma_directory="/path/to/new-references",
                  discriminants=(5,), ranks=(3,), degree=1)
```

`deduplicate_faces=true` applies a recorded performance patch only to the
copied VMH reference source: it avoids expanding different facet tuples with
exactly the same minimum-vector intersection. It supplies no Julia cells or
answers to Magma. The original and patched source hashes are saved in
`source.json`. Omit the option to use the pristine VMH algorithm.

These are expensive optional comparisons, not the default CI suite.
`test/real_quadratic_sl_dimensions.jl` checks determinant-height translation,
odd/even centres, an eight-dimensional collapse, corrupted certificates,
incidence coverage, and agreement of graph links with exhaustive flag links.

Omit the last argument to skip the Magma comparison. Tests and their generated
receipts are independent of the construction; neither expected homology nor
field-specific cell counts are encoded in the construction.

## Detailed construction interface

```@docs
ArithmeticVoronoi.RealQuadraticSL.compute_real_sl_complex
```
