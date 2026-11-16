# Experimental real-quadratic SL₂ complexes

The local `real-quadratic-sl-complex` branch contains a prototype that constructs
an actual three-dimensional contractible SL₂(O)-equivariant cell complex.
It is explicitly loaded from `experimental/RealQuadraticSL`; the normal HAP
provider and registered package entry points are unchanged. It has been exercised
on O for Q(√5) and Q(√3), whose field discriminants are 5 and 12.

## Run a fresh computation

From a checkout with its Julia dependencies instantiated:

```julia
using ArithmeticVoronoi, Oscar
import Oscar: GAP
include("experimental/RealQuadraticSL/RealQuadraticSL.jl")
using .RealQuadraticSL

result = compute_real_sl_complex(5; progress=true);
length.(result.geometric.cells)  # geometric cell orbits, starting at dimension 0
length.(result.chain.cells)      # algebraically reduced permutation modules
K = result.hap;
R = GAP.Globals.FreeGResolution(K, 4);
C = GAP.Globals.TensorWithIntegers(R);
[GAP.gap_to_julia(GAP.Globals.Homology(C,n)) for n in 1:3]
```

The equivalent command-line runner is:

```sh
julia --project=. experimental/RealQuadraticSL/run.jl 5 /absolute/path/to/output
```

Set `ARITHMETIC_VORONOI_HAP_ROOT` to a HAP checkout to select it explicitly.
For constructing the complex, an installed HAP also suffices. No Magma or C++
is called. Enumeration and certification run afresh. The optional output directory
receives link data, exact collapse witnesses and construction timings; these files
are outputs only and are never loaded implicitly. Use `reduce=false` to export
the full geometric complex to HAP instead of its smaller chain-equivalent model.

The prototype currently requires a positive quadratic field discriminant and the
free rank-two lattice over its maximal order. It returns no complex if its chosen
height cannot be certified. This is not a proof that such a complex does not exist
for that field. A contracting homotopy for HAP cup-product calculations is not
provided.

## Construction and certificate

Let W be the bounded Ryshkov cell complex already computed for GL₂(O), and let
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
and GL₂(O) acts by integer translations on height.

The rational b(c) are proposed from approximate determinant heights of invariant
interior forms, with small rational perturbations. Approximation is used only to
choose a candidate. All subsequent comparisons, unit exponents and certificates
are exact.

The construction enumerates the entire ascending and descending link at each
barycentric vertex orbit. It finds elementary simplicial collapses to one point,
then independently replays each witness using the full face set. Greedy collapse
is tried first, followed by finitely many reproducible alternate orders. Failure
to produce a witness aborts construction.

By the Morse lemma and Corollary 2.6 of
[Bestvina–Brady, *Morse theory and finiteness properties of groups*, §2](https://people.math.osu.edu/davis.12/courses/8800-20/Bestvina-Brady.pdf),
these link certificates imply that a regular level has the homotopy type of W.
Thus `h⁻¹(0)` is contractible. SL₂(O) preserves it. The height kernel has
only the additional determinant sign, so SL₂ has finitely many cell orbits
in this level.

The intersection of a mixed-height simplex with the level is a product of two
simplices. Its facets are the mixed-height subflags obtained by deleting one
vertex. Flags have no inversions, hence their stabilizers fix the level cells
pointwise. Exact determinant-one orbit splitting supplies the SL transporters
and stabilizers, including central −I. Boundary signs use the height
coorientation, with zero-cells normalized to positive augmentation.

The implementation verifies equivariant `∂²=0`, augmentation, and invariance of
each boundary under its stabilizer. Optional algebraic cancellations remove pairs
of isomorphic permutation modules; conjugacy of the stabilizers is checked at
every pivot. The full geometric object is retained separately. The reduced
object is a chain model, not a claim that its basis still consists of geometric
product cells. Both can be exported as HAP non-free resolutions.

## Validation and separation from computation

Optional tests are in `test/real_quadratic_sl.jl`; they are not in the default
small test suite. They check certificate replay, reject incomplete witnesses and
changed boundary orientations, and compare H₁ with an independently constructed
SL presentation.

For homology through degree three they separately restrict a GL₂ complex to
`det⁻¹((O×)²) ≅ SL₂(O) × Z`, then check the integral Künneth formula. The same
comparison can use freshly generated Magma/VMH GL₂ complexes, with explicit
conversion of VMH's rational row representation. Magma does not construct the
real-quadratic SL₂ complex in this comparison.

```sh
# Optional: licensed Magma plus a VMH checkout; output must be new.
julia --project=. test/prepare_real_sl_magma.jl /path/to/VMH-TotallyRealNumberFields /path/to/references

# Uses HAP with RestrictNonFreeResolution; no precomputed expected answers.
ARITHMETIC_VORONOI_HAP_ROOT=/path/to/hap julia --project=. test/real_quadratic_sl.jl /path/to/test-output /path/to/references
```

Omit the last argument to skip the Magma comparison. Tests and their generated
receipts are independent of the construction; neither expected homology nor
field-specific cell counts are encoded in the prototype.
