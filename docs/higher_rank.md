# Higher-rank GL presentations and HAP complexes

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

# A full presentation and constructive membership for GL(3,Z[i]).
R = run_gl_pipeline(voronoi_request(-4; rank=3);
    directory="results/my_GL3_gaussian");
S = R.simplified[:GL];
F = free_group(length(S.generator_matrices));
G, toG = quo(F, [F(r) for r in S.relators]);

P = R.enumeration.problem;
g = identity_matrix(P.field, 3); g[1,3] = P.generator;
c = matrix_to_word(R, g);              # GL is the default for this runner
verify_word_certificate(R, g, c)      # true, exact geometric path replay
R.run["presentation"]                # sizes and abelianization
```

Use `run_gl_pipeline(classical_problem(; rank=4); directory="results/my_GL4_Z")`
for GL₄(Z), or `run_gl_pipeline(5; rank=3, directory="results/my_GL3_D5")`
for GL₃ over the maximal order of Q(√5). The rank is not capped, but the cost
of enumerating forms, stabilizers and ridges grows rapidly.

The runner writes the enumeration, complete attaching-loop certificates,
raw and simplified presentations, Tietze maps and stage timings. Compatible
saved loops are exactly replayed on reruns. It uses the native enumerator and
imposes no timeout or memory limit. `simplify=false` skips Tietze simplification;
then use `matrix_to_word(R, g; simplified=false)`.

In higher rank it walks representatives under each full vertex stabilizer.
The saved coverage ledger includes a checked matrix and stabilizer word
identifying every omitted ridge with its representative. It computes a
cone's extreme rays once and restricts them to its ridge faces. These
reductions preserve complete coverage without repeating equivalent work.

Run quadratic request records in one Julia process:

```sh
julia --project=. scripts/run_gl.jl --requests=examples/higher_rank_gl_requests.json --output=results/higher_gl
```

For higher homology, construct the full equivariant complex and pass it to HAP:

```julia
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

`run_gl_pipeline` constructs GLₙ presentations using the 2-skeleton in every
supported rank. The combined GL/SL/PSL `run_field_pipeline` remains rank two,
as does the optional C++ request adapter. They reject higher-rank requests
explicitly. Rational and
real-quadratic SL complexes, projective quotients, nonmaximal orders, nonfree lattices and
fields of degree greater than two remain outside the provider's scope.
The record adapter takes quadratic fields; use `classical_problem` for Q.

Supporting an input rank does not imply practical computation at arbitrary
rank: the number of forms, faces and transporter candidates grows rapidly.
HAP's finite-stabilizer resolutions and integral homology computations also
have costs beyond constructing the cell complex.

### What HAP still needs

A GL presentation is sufficient for group operations and H₁, but does not
replace the equivariant complex used by HAP's higher-homology path. That path
needs cells in the relevant dimensions, their full finite stabilizers,
orientation actions and boundary transporters. The native provider already
supplies these in higher rank; no further HAP dispatcher change is needed.
The provider does not supply a contracting homotopy for the resulting free
resolution. Features requiring one, including some cup-product routines, need
additional work. Generating GL presentations only requires the 2-skeleton
and does not require computing H₂.

## Regression tests and independent comparisons

The GL presentation regressions construct GL₂(Z), GL₃(Z) and GL₃(Z[i])
from scratch, check all matrix relations, replay words for elementary,
permutation and unit matrices, and reject damaged word certificates. For
the rank-three cases they also compare with presentations extracted by
HAP's `PresentationOfResolution` and with integral H₁. They are included in
`test/hap.jl` and can be run separately with `test/gl_presentations.jl`.

To compare GL presentations with fresh Magma computations:

```sh
ARITHMETIC_VORONOI_HAP_ROOT=/path/to/hap julia --startup-file=no --project=. scripts/compare_gl_presentations.jl /path/to/VMH-TotallyRealNumberFields /path/to/VMH-ImaginaryQuadraticNumberFields results/gl_validation
```

The comparison checks a bijection of the perfect forms by exact isometries
and uses those isometries and the native coherent vertex lifts to align
coordinates. The imaginary VMH matrices are converted by
`h^-1 * transpose(g^-1) * h`; the real VMH code already uses column matrices
and uses `h^-1 * g * h`. These are group automorphisms, with the field
embedding and h saved alongside the maps.
Exact matrix search finds short generator maps, and Julia's geometric word
certificates are retained too. A second Magma phase
checks that both maps respect the abstract relators and are inverse on the
abstract generators; its completion is recorded separately from matrix
checks. It uses derived relation rewrite rules: ordinary Magma fp-word
equality only does free reduction. An inconclusive reduction expands the
completion search and tries other generator orderings; it is never reported
as a mathematical counterexample. Only positive reductions to the identity
are accepted as proofs.
The reference source trees are copied per case and their algorithms are
unchanged. The tests below additionally compare
the full equivariant complexes and homology.

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

The earlier [complex validation receipt](../results/higher_rank_validation.json)
records the tested versions, source hashes and completed comparisons. That
run passed 1,329 native assertions and 239 HAP assertions with bundled
HAP 1.78, or 248 with the pinned upstream HAP 1.80 checkout. Its rank-four
rational and rank-three real-quadratic checks compared perfect forms only;
presentation validation is recorded separately below.

The completed complex and homology comparisons in that receipt are:

| Group | Perfect forms | H₁ | H₂ |
| --- | ---: | --- | --- |
| GL₃(Z) | 1 | C₂ | C₂ × C₂ |
| GL₂(O₅) | 2 | Z × C₂ | C₂⁴ |
| GL₂(O₋₇) | 1 | C₂³ | C₂ × C₁₂ |
| SL₂(O₋₇) | 1 GL orbit | Z × C₄ | C₁₂ |
| GL₃(Z[i]) | 1 | C₄ | extended comparison incomplete |
| SL₃(Z[i]) | 1 GL orbit | trivial | extended comparison incomplete |
| GL₃(O₅) | 7 | outside that run | outside that run |
| GL₄(Z) | 2 | outside that run | outside that run |

Here O_D is the maximal order of the quadratic field of discriminant D, and
C_m is cyclic of order m. The Gaussian GL₃/SL₃ cell counts, stabilizer orders
and H₁ also agree with the fresh Magma complexes. The receipt records the
incomplete extended H₂ comparison separately.

### GL presentation validation

The [presentation validation receipt](../results/gl_presentation_validation.json)
records four fresh comparisons with Magma 2.30-1. Each checks exact perfect-form
isometries and generator matrices, then proves reciprocal homomorphisms
between the abstract presentations. Thus these comparisons establish more
than agreement of abelianizations or matrix relations.

| Group | Perfect forms | Simplified generators / relators | H₁ | Magma presentation isomorphism |
| --- | ---: | ---: | --- | --- |
| GL₃(Z) | 1 | 3 / 10 | C₂ | verified |
| GL₄(Z) | 2 | 3 / 63 | C₂ | verified |
| GL₃(Z[i]) | 1 | 3 / 15 | C₄ | verified |
| GL₂(O₋₇) | 1 | 4 / 9 | C₂³ | verified |
| GL₃(O₅) | 7 | 3 / 207 | Z × C₂ | reference calculation incomplete |

HAP's independently assembled presentations and H₁ were also checked for
GL₃(Z), GL₃(Z[i]) and the GL₂(O₋₇) regression. The new presentation tests
pass 122 assertions with each of HAP 1.78 and 1.80; the native regression
tests pass 1,329 assertions, for 1,573 assertions in the recorded runs.
These runs use Julia 1.13.2 and OSCAR 1.8.2.

The GL₃(O₅) Julia calculation is complete: its 11,426 local ridges reduce to
391 stabilizer orbits and 782 oriented attaching loops, with exact coverage
witnesses for all ridges. The pristine Magma reference independently found
seven perfect forms but was stopped manually during its exhaustive
attaching-loop calculation. The earlier receipt verifies the perfect-form
bijection; no independent presentation isomorphism is claimed for this case.
The comparison runner retains this case so it can be completed separately.

To run only that supplementary comparison (it can take much longer than
the native calculation):

```sh
ARITHMETIC_VORONOI_HAP_ROOT=/path/to/hap julia --startup-file=no --project=. -e 'include("scripts/compare_gl_presentations.jl"); compare_gl_presentations(ARGS; cases=((D=5,rank=3,hap=false),))' /path/to/VMH-TotallyRealNumberFields /path/to/VMH-ImaginaryQuadraticNumberFields results/gl_real_validation
```
