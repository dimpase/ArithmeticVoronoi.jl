# Certified real-quadratic SLₙ complexes

The determinant-height construction accepts every rank `n ≥ 2` over a real
quadratic maximal order O, on the free lattice Oⁿ. It constructs an
SLₙ(O)-equivariant regular level of the GLₙ Voronoi complex and certifies its
contractibility by exact ascending and descending link collapses. The full
level has dimension `n²−1`. Every run must pass certification; an accepted rank
is not a promise that enumeration or certification will finish quickly.

Rank two was introduced in ArithmeticVoronoi 0.1.1. Higher ranks require the
subsequent development version. Validation covers integral homology through
degree three for SL₂ over Q(√5) and Q(√3), and through degree one for SL₃ over
Q(√5). Higher ranks use the same construction and certification, but have not
yet been run through the full pipeline; see the validation details below.

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
are outputs only and are never loaded implicitly. Use `reduce=false` to retain
the geometric cells. The default algebraic reduction preserves the claimed
resolution degree; the additional truncation for proper partial models is
described below.

For rank three, the corresponding calls are:

```julia
K = hap_voronoi_complex(5; rank=3, group=:SL, resolution_length=2, progress=true);
R = GAP.Globals.FreeGResolution(K, 2);
C = GAP.Globals.TensorWithIntegers(R);
GAP.gap_to_julia(GAP.Globals.Homology(C,1))  # []: trivial H₁
# Or, after enable_hap_julia_backend(; resolution_length=2):
K = GAP.evalstr("VoronoiComplexSL(3,5 : backend := \"julia\")");
```

`resolution_length=m` computes only the cells through dimension `m` when the
full level is larger. It still certifies **all** ascending and descending links
of the full level. The returned HAP object supports a free resolution through
that degree, hence homology through degree `m−1`. Its
`arithmeticVoronoi.computedCellDimension` and `certifiedCellDimension` metadata
record the distinction. Asking it for uncomputed higher cells raises an error.
A later call with a larger resolution length computes a larger model afresh.
The default section model retains a certified partial level subcomplex before
taking its skeleton; the connectivity argument is described below.

The detailed `compute_real_sl_complex` interface constructs the full level by
default. Set `max_cell_dimension=m` explicitly to request a partial model; then
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
are exact. An exact sweep over the edge-flag orbits then translates all heights
by a common rational constant to minimise the number of SL vertex orbits in the
regular level. This preserves every ascending/descending comparison. The chosen
offset and before/after vertex counts are recorded in the construction metadata.
Flag enumeration skips branches whose entire height interval misses the level.

The construction enumerates the entire incidence poset at each barycentric
vertex orbit. Comparable cells form edges; chains form simplices, so each link
is the flag complex of its graph. Its ascending and descending subcomplexes
are induced flag complexes. This avoids expanding every maximal flag.

The certificate records dominated-vertex and dominated-edge removals. Each
removal has a cone link, checked by exact neighbour-set inclusion, and expands
into elementary simplicial collapses. An independent replay checks every cone
condition and requires exactly one vertex to remain. See
[Boissonnat–Pritam, *Edge Collapse and Persistence of Flag Complexes*](https://doi.org/10.4230/LIPIcs.SoCG.2020.19).
The search can remove entire orbits under a cell stabiliser. Cell dimensions
are invariant colours: a removed vertex has a different colour from its apex,
and a removed edge uses two colours different from its apex. Thus members of
a removal orbit cannot delete the cone data needed by another member. Every
orbit removal is recorded once. An independent replay verifies that the proposed
permutations are colour-preserving graph automorphisms, checks the representative
cone, and reconstructs its complete removal orbit. Proper dimension colours ensure
that orbit members cannot remove one another's apices or cone edges. The certificate
therefore describes the same elementary collapses without listing every equivalent
step. Tests also expand these compact witnesses and replay individual removals.
Saved certificates carry a format version and distinguish orbit steps from ordinary
graph steps. No numerical topology test is substituted for the certificate.
When cone removals stop, the search also tests vertex links for a complete graph
collapse. A successful witness collapses the vertex's star onto its link: cone
each inner collapse pair with the removed vertex, then collapse the final edge.
The replay checks the entire inner witness before removing the vertex orbit.
Tests expand this relative collapse into elementary simplicial pairs as well.
Failure to produce a witness aborts construction.

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
of isomorphic induced orientation modules; stabilizer conjugacy and compatibility
of their orientation characters are checked at every pivot. The full geometric object is retained separately. The reduced
object is a chain model, not a claim that its basis still consists of geometric
product cells. Both can be exported as HAP non-free resolutions.

## Certified cell sections

The default `cell_model=:sections` avoids many barycentric cells while keeping
an exact certificate. The detailed interface also accepts `cell_model=:flags`
to retain every intersecting flag cell. For example:

```julia
result = ArithmeticVoronoi.RealQuadraticSL.compute_real_sl_complex(5;
    rank=3, cell_model=:sections, resolution_length=2,
    max_cell_dimension=2, progress=true);
K = result.hap;
```

For an original d-dimensional cell, radial projection identifies a level
component with a regular neighborhood of its opposite-sign boundary core.
The core is a full induced subcomplex of the subdivided boundary sphere.
If it collapses to one vertex, its neighborhood is a PL ball, so the component
can be retained as one cell. This uses the regular-neighborhood criteria in
[Buoncristiano, Part I §4.5](https://www.maths.gla.ac.uk/~mpowell/sandro.pdf).
Every proposed collapse is independently replayed. A failed ball proof triggers
subdivision of that section and its boundary; it never permits an unproved merge.

A bounded-degree calculation can omit some higher-dimensional sections as well.
Suppose a core collapses onto a residual flag complex of dimension r. The
implementation verifies the complete partial collapse and computes the residual
clique dimension exactly. Its regular neighborhood N is still a neighborhood
of the remaining core. General position and the mapping-cylinder description
of regular neighborhoods imply that `(N, ∂N)` is `(d−r−2)`-connected; see
[Bryant, Theorems 3.5 and 3.16 and §4](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/pltop.pdf)
and the dimension argument in
[Ferry, Lemma 10.12](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/ferrygt.pdf).
This is an application of those results to the level sections.

For a resolution through m, the code chooses an original-cell dimension D
so that every omitted section has `d−r−2 ≥ m`, and retains every section of
original cells through D. Adding the omitted pieces cannot change homotopy
below m. Thus the retained subcomplex is `(m−1)`-connected, and its m-skeleton
provides the exact augmented cellular chains needed through resolution degree m.
It need not be contractible. The full level's contractibility certificate is
still required, and requests beyond the computed resolution degree are rejected.
The HAP metadata records `levelModel`, `retainedOriginalDimension`,
`computedCellDimension` and the full level's `certifiedCellDimension`.

Ball sections may have orientation-reversing stabilizers. Their actual even
stabilizers are passed to HAP. The boundary validator and chain cancellations
use the orientation characters, including signs from reducing transporter
cosets. Subdivided sections use the original flag-cell orientations. Geometric
cells and the reduced chain model remain separate objects.

For a proper partial model through degree m, the reduced chain model also omits
top generators whose boundary is zero or exactly repeats that of a retained
generator with the same orientation module. Each omission records its source
and, for a repeat, the retained source in `omitted_top_generators`; these indices
refer to `chain.original.cells[end]`. This preserves
the image of the top boundary map exactly, and therefore preserves homology below
m and the claimed resolution degree. It can change the unclaimed top homology;
it is not a chain homotopy equivalence in all degrees. The geometric complex
retains every cell. Full-dimensional models are left intact by this step.
HAP metadata reports the number of `omittedTopGenerators`.

## Validation and separation from computation

The HAP integration CI tests `hap_voronoi_complex` and the registered GAP
provider on discriminants 5 and 12, checking integral homology through degree
three. These tests are in `test/hap_real_quadratic_sl.jl`.
`test/real_quadratic_sl_sections.jl` additionally compares coarse, mixed and
fully subdivided cell models, checks orientation-aware cancellation, and rejects
requests beyond the certified partial degree. Both run in HAP CI.
`test/hap_orientation_characters.jl` checks the adapter's shortcut for proved
trivial orientation characters against the general coset formula, including
transporters outside the stabilizer and an orientation-reversing example.

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
                     cases=((D=5,rank=3),), deduplicate_faces=true, max_cell_dimension=2,
                     memoize_cell_data=true, linear_tangent_spaces=true, incremental_boundaries=true)
include("test/real_quadratic_sl.jl")
run_real_sl_tests("/path/to/new-test-output"; magma_directory="/path/to/new-references",
                  discriminants=(5,), ranks=(3,), degree=1)
```

`deduplicate_faces=true` applies a recorded performance patch only to the
copied VMH reference source: it avoids expanding different facet tuples with
exactly the same minimum-vector intersection. It supplies no Julia cells or
answers to Magma. The original and patched source hashes are saved in
`source.json`. `memoize_cell_data=true` also avoids recomputing a cell's Magma
stabilizer and affine tangent basis for every orientation query. The calculations
and returned values are unchanged; cached values are produced by Magma within
that run. The optional `linear_tangent_spaces=true` computes the same affine
tangent space as the common kernel of the minimum-vector evaluation equations, using VMH's
`FindAllPerp`. This avoids enumerating boundary embeddings for orientation queries;
the subsequent boundary computation uses VMH's original algorithm unless
`incremental_boundaries=true` is also requested. Its
orientation signs have been checked against the original calculation for every
stabilizer element of the full rank-two reference complexes at discriminants 5
and 12 (318 comparisons). All performance options are opt-in. Source changes, including
the dimension bounds below, are recorded with original and modified hashes.
The orientation comparison is reproducible with an original prepared reference:

```sh
julia --project=. test/compare_real_sl_magma_tangents.jl /path/to/original-references/magma_D5 /path/to/new-orientation-check
```

`incremental_boundaries=true` searches the same ordered independent images of
VMH's chosen minimum-vector basis, testing each linear relation as soon as its
last nonzero coordinate is assigned. A rejected prefix cannot extend to an
embedding. The integral-unit test and stabilizer-coset test are unchanged, and
all vectors still come from Magma. Exact comparisons with the original search
passed for all 331 adjacent-class pairs in the two rank-two references, and 16
independently constructed positive-form pairs in rank three. To reproduce them:

```sh
julia --project=. test/compare_real_sl_magma_boundaries.jl /path/to/original-references/magma_D5 /path/to/new-boundary-check
# A rank-three prepared reference can test small forms without a full enumeration:
julia --project=. test/compare_real_sl_magma_boundaries.jl /path/to/original-references/magma_rank3_D5 /path/to/new-rank3-check --configurations-only
```

The optional Magma runners write their completion marker inside the procedure
that performs the calculation, so a failed assertion or runtime error cannot
be followed by a success marker.

`max_cell_dimension=m` explicitly limits the Magma reference to its m-skeleton.
The comparison requires `m ≥ degree+1`; it never treats the skeleton as a
contractible full complex. This skips unneeded higher-dimensional classes and
is recorded with the reference source hashes. Omit it to compute all dimensions.

These are expensive optional comparisons, not the default CI suite.
`test/real_quadratic_sl_dimensions.jl` checks determinant-height translation,
odd/even centres, an eight-dimensional collapse, corrupted certificates,
incidence coverage, and agreement of graph links with exhaustive flag links.

Omit the last argument to skip the Magma comparison. Tests and their generated
receipts are independent of the construction; neither expected homology nor
field-specific cell counts are encoded in the construction.

## Rank-three validation

For SL₃ over the maximal order of Q(√5), the certified partial level through
resolution degree two gives **H₁(SL₃(O), Z) = 0**. The calculation independently
replayed all 3,878 ascending and descending link certificates and all 3,718
section certificates. Its connectivity bounds required original GL₃ cells
through dimension five; that cutoff was computed from the certificates.

| Model | 0-cells or generators | 1-cells or generators | 2-cells or generators |
|---|---:|---:|---:|
| Geometric partial level | 1,142 | 8,432 | 20,399 |
| After equivariant cancellations | 11 | 35 | 13,133 |
| After preserving only the required top boundary image | 11 | 35 | 1,514 |
| HAP free resolution | 11 | 57 | 1,594 |

The HAP result agrees with a separately constructed SL₃ presentation and with
the determinant restriction of the native GL₃ complex, whose first homology is
Z. An independent fresh Magma/VMH GL₃ computation also agrees: its first three
cell-orbit counts are `[7,34,140]`, all corresponding stabilizer-order multisets
match, and HAP gives H₁ = Z after restriction to SL₃ × Z. This comparison used
VMH's original GAP reader; the optional Magma source performance changes above
were recorded with source hashes and supplied no Julia cells or expected answers.
A separate Nemo check of the exported integer boundary matrices verified
that their product is zero, H₀ is Z and H₁ is zero using Smith normal forms.
Both the general orientation-character formula and the optimized adapter gave
the same result. These checks establish H₁ for this example; they do not claim
higher homology in rank three or a completed example in rank four.

On the validation host, equivariant cancellation took 81 seconds. HAP's free
resolution took 29 seconds, evaluating its integer boundary matrices took
31 seconds, and the final homology computation took 0.25 seconds. The independent
Magma reference took 202 CPU seconds to enumerate its perfect forms and 4,579
CPU seconds to construct its two-skeleton. Julia used OSCAR 1.8.2 and HAP 1.80;
the reference used Magma 2.30-1. Enumeration,
link and section certification are additional costs in the SL calculation;
these are phase timings,
not a complete first-call benchmark. Development validation explicitly reused
SHA-verified intermediate outputs between phases. Every geometric input and
certificate was computed from scratch, and the package computation has no
implicit snapshot or expected-answer loader.

## Detailed construction interface

```@docs
ArithmeticVoronoi.RealQuadraticSL.compute_real_sl_complex
```
