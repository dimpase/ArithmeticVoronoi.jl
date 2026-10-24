# Higher-rank native enumeration and HAP complexes

The native enumerator and HAP provider accept ranks `n >= 2`. Quadratic
requests describe GLₙ over a real or imaginary quadratic maximal order, or
SLₙ over an imaginary quadratic maximal order, on the free lattice Oⁿ.
Higher-rank support uses the same general record interface merged in
[HAP PR #167](https://github.com/gap-packages/hap/pull/167); HAP needs no new
dispatcher to accept a larger rank.

## Compute in Julia

From the checkout root, start `julia --startup-file=no --project=.`:

```julia
using ArithmeticVoronoi, Oscar
import Oscar: GAP

# GL(3,Z[i]); -4 is the field discriminant.
request = voronoi_request(-4; rank=3);
K = hap_voronoi_complex(request; resolution_length=2, progress=true);
R = GAP.Globals.FreeGResolution(K, 2);
C = GAP.Globals.TensorWithIntegers(R);
sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(C, 1)))) # [4]
```

Use `group=:SL` in the request to select SL₃. To enumerate forms first and
inspect the cells:

```julia
E = enumerate_perfect_forms(request; progress=true);
V = voronoi_complex(E; group=:GL, progress=true);
length(E.representatives)  # 1
length.(V.cells)           # [1,1,3,5,4,3,2]
validate_voronoi_complex(V)
```

The lower-level model also supports rational forms:

```julia
E = enumerate_perfect_forms(classical_problem(; rank=3); backend=:native);
V = voronoi_complex(E);
length.(V.cells) # [1,1,2,1], for GL(3,Z)
```

To run the Gaussian GL₃ and SL₃ examples in one Julia process:

```sh
julia --project=. scripts/run_hap.jl --requests=examples/higher_rank_hap_requests.json --degree=1 --output=results/rank3_homology
```

The batch writes `GL3_D-4.json`, `SL3_D-4.json` and `summary.json`. Rank-two
filenames retain their old names. `--resume` compares the entire request,
including rank and homology degree.

## Scope

The form-space basis, exact trace coordinates, finite stabilizers, neighbour
traversal, spanning tests, cell transporters, orientations and equivariant
boundaries use the requested rank. All native complexes are checked for
`d²=0` before passing to HAP. SL uses the full determinant kernel, with no
quotient by scalar units. In odd rank, −I is excluded from SL; other central
roots of unity with determinant one must still be retained.

The GL/SL/PSL **presentation pipeline remains rank two**, as does the optional
C++ request adapter. They reject higher-rank requests explicitly. Rational and
real-quadratic SL complexes, projective quotients, nonmaximal orders, nonfree lattices and
fields of degree greater than two remain outside the provider's scope.
The record adapter takes quadratic fields; use `classical_problem` for Q.

Supporting an input rank does not imply practical computation at arbitrary
rank: the number of forms, faces and transporter candidates grows rapidly.
HAP's finite-stabilizer resolutions and integral homology computations also
have costs beyond constructing the cell complex.

## Regression tests and independent comparisons

Run the native HAP integration tests, including the existing rank-two cases:

```sh
julia --startup-file=no --project=. test/hap.jl
```

To test a specific HAP checkout, select it before HAP is loaded:

```sh
ARITHMETIC_VORONOI_HAP_ROOT=/absolute/path/to/hap julia --startup-file=no --project=. test/hap.jl
```

This mode checks HAP's own dispatcher before enabling the Julia adapter and
verifies that registration does not replace its entry points. CI is configured
to run both bundled HAP and a pinned upstream checkout containing the merge.

Higher-rank tests cover trace models through rank four, the two classical
rank-four perfect forms, rational rank-three H₁/H₂, and Gaussian rank-three
GL/SL H₁. The Gaussian GL result is checked against HAP's stored
`ResolutionGL3QuadraticIntegers(-1,2)`. A separate Eisenstein rank-three
stabilizer test checks that SL retains its centre of order three and excludes −I.

For fresh Magma comparisons, obtain these independent implementations:

- [VMH-TotallyRealNumberFields](https://github.com/schoennenbeck/VMH-TotallyRealNumberFields),
  revision `befcfeda041781b7de922a16a37388decee2c642`.
- [VMH-ImaginaryQuadraticNumberFields](https://github.com/schoennenbeck/VMH-ImaginaryQuadraticNumberFields),
  revision `b5b79a8103c648c36e90b8bd95880782a3ff0d2e`.

With a working licensed `magma` on PATH:

```sh
julia --startup-file=no --project=. scripts/compare_magma_higher_rank.jl /path/to/VMH-TotallyRealNumberFields /path/to/VMH-ImaginaryQuadraticNumberFields results/higher_rank_validation
```

The runner copies each checkout for each case and changes only `BasicData.m`.
It checks fresh Magma perfect forms by exact K-linear isometries, then compares
cell counts, stabilizer orders and homology of the independent complexes
through HAP: H₁/H₂ for the rational rank-three and quadratic rank-two cases,
and H₁ for Gaussian rank three. GL H₁ is also compared with Magma's
independently computed presentation.
The cases are rational and Gaussian rank three, real D=5 and imaginary D=−7
rank-two regressions, and perfect-form enumeration for real D=5 in rank three
and for Q in rank four.
Magma inputs, logs, complexes and a JSON comparison receipt remain in the
output directory. No Magma installation is needed for the regular tests.

Appending `--extended` also requests Gaussian rank-three H₂. This independent
comparison is still incomplete: the validation attempts were stopped manually
during HAP's homology computation. It is not part of the passing regression
claim, and the scripts impose no timeout or memory limit.

Completed independent enumeration comparisons include seven perfect forms
over Q(√5) in rank three and the two rational rank-four forms, with exact
isometries giving a bijection between the Julia and Magma representatives.
These enumeration checks do not assert that the corresponding higher-rank
group presentations or homology have been computed.

The [local validation receipt](../results/higher_rank_validation.json) records
the tested versions, source hashes and completed comparisons. The native
regression suite passes 1,329 assertions; the HAP suite passes 239 with bundled
HAP 1.78 and 248 with the pinned upstream HAP 1.80 checkout.

The completed fresh Magma comparisons include:

| Group | Perfect forms | H₁ | H₂ |
| --- | ---: | --- | --- |
| GL₃(Z) | 1 | C₂ | C₂ × C₂ |
| GL₂(O₅) | 2 | Z × C₂ | C₂⁴ |
| GL₂(O₋₇) | 1 | C₂³ | C₂ × C₁₂ |
| SL₂(O₋₇) | 1 GL orbit | Z × C₄ | C₁₂ |
| GL₃(Z[i]) | 1 | C₄ | extended comparison incomplete |
| SL₃(Z[i]) | 1 GL orbit | trivial | extended comparison incomplete |
| GL₃(O₅) | 7 | not computed | not computed |
| GL₄(Z) | 2 | not computed | not computed |

Here O_D is the maximal order of the quadratic field of discriminant D, and
C_m is cyclic of order m. The Gaussian GL₃/SL₃ cell counts, stabilizer orders
and H₁ also agree with the fresh Magma complexes. The receipt records the
incomplete extended H₂ comparison separately.
