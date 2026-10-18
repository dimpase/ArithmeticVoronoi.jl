# Constructive matrix words: implementation and validation scope

The existing cone-descent routine was inspected and retained. It is now
exposed through a public GL/SL/PSL wrapper with exact path replay. The supported
fields are discriminants 5, 8, 12, 13 and 40, on the free maximal-order lattice
O^2. Within each completed case the input may be an arbitrary valid matrix;
it is not restricted to words in the final generators or to scalar units.
No claim of support for other field models is made by the current constructor.

## Public interface

```julia
using ArithmeticVoronoi, Oscar
R = run_field_pipeline(12)
P = R.enumeration.problem
T = matrix(P.field, [1 P.generator; 0 1])
c = matrix_to_word(R, T; target=:SL, cache_path="tau_w.json")
c["word"]
verify_word_certificate(R, T, c)
```

`target` is `:GL`, `:SL` or `:PSL`; its default is `:SL`. By default the word
uses the corresponding saved simplified presentation. `simplified=false`
selects the raw alphabet; `certificate=false` returns only its signed-index
word. The word product is evaluated in written matrix-product order. It is
neither a shortest word nor a canonical normal form.

A two-by-two matrix is interpreted over the stated K. An encoded four-by-four
integral matrix is also accepted only if it commutes with multiplication by
w and decodes to the same K-linear action. Integral entries in O and unit
determinant are checked before GL reduction. SL/PSL require determinant one.
For PSL the input is an SL lift: the returned lift may differ by -I, with
its exact sign in `SL_lift_sign`. This does not accept arbitrary PGL matrices.

The original method `matrix_to_word(C::PresentationCertificate,g)` remains a
raw GL interface. The new wrapper uses the existing Schreier rewrite to move
determinant-one words through the fundamental-scalar quotient into SL, then
uses the saved Tietze maps. The determinant-one scalar correction is unique,
so the SL result equals the input matrix including its sign.

## Why the reducer is general and terminates

The root interior B is a positive sum of every root-cone evaluation ray.
The target Q=g^-1 B g^-t lies in the interior of the unique target full cone.
The code examines every exact facet of the current representative and crosses
an actually separating facet; it is not a breadth-first search in a fixed
small set of actual cones. It updates frames by b*s_src*a and prepends the
corresponding coherent words.

At each crossing the exact objective Tr_K/Q(tr(F_current Q)) decreases.
Let N_F clear the field-coordinate denominators of all finitely many
representative forms, and N_Q clear those of Q. Because all frames are
O-integral and the field trace of an algebraic integer is integral, every
objective lies in (1/(N_F N_Q)) Z. The objective is positive, hence there are
only finitely many strictly decreasing values below its initial value.
`objective_denominator` records this bound, and every step checks it.
Positivity and decrease without this fixed denominator would not suffice.

On termination F_root^a=F_root^g. The finite root-stabilizer correction is
h=a*g^-1 and the final word is h^-1*a. The certificate retains its word in
the root finite-group alphabet; -I is never replaced by the identity of an
unsigned permutation action. The source correction S_i and target reverse
correction R_j are distinguished in the proof and convention regression.

## Certificates, failure and resource scope

A complete certificate records its field/fan source hash, alphabet fingerprint,
exact input, target group, final word and lift sign. Its raw GL path records
source/target cones, source and next frames, both crossing labels, facet IDs,
trace objectives and denominator bound, accumulated words and finite residual.
Replay recomputes each crossing without searching and checks the final group
translation. It rejects changed input, alphabet, source, frames or corrections.

`cache_path` atomically stores a complete certificate. Loading it replays and
verifies it; it is not trusted merely because its matrix word evaluates
correctly. A step/time budget exhaustion creates an incomplete `.partial`
checkpoint and raises `WordReductionError(:resource_exhausted,...)`.
Incomplete paths are retained diagnostically; a new search currently restarts
at the root rather than resuming a partly traversed path. Complete caches
are reusable. `max_steps` limits reduction iterations, including the terminal
cone check. `timeout` is checked between exact crossings; it is not a hard
preemption of a single arithmetic operation or of Julia startup/compilation.

Invalid membership/budgets raise `:invalid_input`; an unavailable simplified
alphabet raises `:unavailable_alphabet`; a corrupted cache raises
`:invalid_certificate`. Other reduction failures remain explicit failures and
never produce a complete word. There is no presentation or word claim for
an unfinished case.

## Validation evidence

`test/fields_words.jl` exercises independently specified identity, -I, finite
stabilizer elements, eI/e^-1 I, diag(e,e^-1), diag(e,1), upper/lower elementary
translations, and deterministic products with positive and negative powers.
It checks fresh paths and cache replay, all three target groups, raw and
simplified alphabets, K-linear encoded inputs, a path leaving the root cone,
central signs, invalid determinants/denominators/semilinearity, explicit
exhaustion, and corrupt frames and residual corrections. It is run together
with all original regressions. The authoritative completion/status evidence
is `results/fields_words_tests.json` and the fresh test transcript
`results/fields_delivery_validation.log`.

The D8 hand fixture and four-field Kirchheimer-Wolfart maps give additional
independently specified matrix correspondences. They do not replace the
termination argument or geometric completeness proof. The D40 nonprincipal
ideal regression and full pipeline ensure the general reducer has no hidden
class-number-one or norm-Euclidean division assumption. The reference-only
row-reduction helper is excluded from the D40 computation.
