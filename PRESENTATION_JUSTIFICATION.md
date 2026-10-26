# Justification of the arithmetic Voronoi presentations

The original computations in `results/two_cells/D5/` and `D8/`, and the
extended pipeline in `results/fields/D*/`, concern the full K-linear groups
GL₂(O), SL₂(O), and the separately constructed PSL₂(O), for maximal orders
of real quadratic fields with discriminants 5, 8, 12, 13 and 40.
Only records whose stage status is complete claim a presentation.
GL and SL retain −I in these rank-two examples; PSL is constructed separately.
No Galois semilinear extension is included. This document explains
the completeness argument; testing matrix relations alone would not establish
it. The certificates are exact computational certificates, not proof-assistant
formalizations.

The argument below is written for real quadratic forms and the original
five-field certificates. The current native implementation also covers the
[16-field real survey](docs/native_survey.md) and
[imaginary quadratic Hermitian forms](docs/imaginary_quadratic.md).
The GL, SL and PSL constructions are now rank-independent; the changes of dimension and
the supported arithmetic scope are specified below.

The reference is Coulangeon, Nebe, Braun and Schönnenbeck,
[*Computing in arithmetic groups with Voronoï's algorithm*](https://arxiv.org/abs/1407.6234),
version 1: Lemma 2.3 and Theorem 2.4 (perfect cones), Section 3 and
Theorem 3.5 (the contractible well-rounded complex and duality), Theorem 4.1
(the group-action presentation theorem), and Section 7.3 (ridge cycles).
Below we specify how the hypotheses and path conventions are realized here.

## General GL rank

`run_gl_pipeline` applies the same construction to Mₙ(K) acting on Oⁿ,
for every integer n ≥ 2, over Q or a real or imaginary quadratic field.
The order is maximal and the lattice is free. `run_field_pipeline` applies
the determinant-kernel and projective-quotient constructions below to the
resulting GL presentation in the same rank.

Write m for the dimension of the rational form space: n(n+1)/2 over Q,
n(n+1) for a real quadratic field, and n² for an imaginary quadratic field.
Replace the dimensions six, five and four below by m, m−1 and m−2.
The trace pairing is formed in the entire m-dimensional space. Supporting
vectors of a positive ridge must span Kⁿ; exact column-rank computation
records n independent vectors. Their sum of outer products is positive
definite at every embedding (using conjugate transpose in the Hermitian
case), and exact positivity is checked. Nonspanning ridges are boundary
faces and are discarded: they do not define attaching 2-cells of the
well-rounded complex. In particular, GL₂(Z) has a graph and no such 2-cells.

The transverse ridge link is still a circle. Ordered crossings, reverse
replay, complete ridge coverage, full edge stabilizers, coherent tree lifts
and the group-action theorem are unchanged. Only the 2-skeleton is needed
to present the group. All matrices and constructive-word denominator
checks use the requested field and rank.

The higher-rank runner walks both orientations of one ridge in each orbit
of the full finite vertex stabilizer. For every omitted local ridge, the
coverage ledger records a stabilizer matrix, its word in the original
stabilizer generators, and the chosen representative. Validation replays
the word, tests membership in the full vertex stabilizer, and checks its
exact action on the complete supporting-vector set. The ledger must cover
every local ridge, and each chosen representative must have both attaching
cycles. These are translates of the same attaching disk under an already
included vertex stabilizer, so the selected cycles still meet every global
2-cell orbit. No global 2-cell orbit count is asserted. The historical
rank-two runner continues to walk every local ridge individually.

Extreme rays of every ridge are obtained by restricting the parent's exact
extreme-ray set to the face. This is the face-of-a-polyhedral-cone property;
the independent face dimension and positivity checks remain in place.

Higher-rank finite stabilizers store every element, its normal-form word,
and right multiplication by every chosen generator. The relations
`w_x g = w_(xg)` still give a complete finite-group presentation: a full
order-by-order multiplication table adds no relations needed for this
argument. The serialized `generator_products` table replaces the full
`multiplication_table` in higher-rank certificates.

## SL and PSL in arbitrary rank

For Q and imaginary quadratic K, the unit group O× is finite cyclic of order
m. The determinant maps GLₙ(O) onto O×. Reidemeister–Schreier rewriting of
its kernel, with index m, gives SLₙ(O) directly. No scalar unit is killed
at this stage.

For real quadratic K, write O× = {±1} × ⟨ε⟩, with ε a verified fundamental
unit. Constructive reduction gives a GL word for εI. After imposing εI=I,
determinant induces a surjection

    GLₙ(O)/⟨εI⟩ → O×/⟨εⁿ⟩ ≅ C₂ × Cₙ.

Its kernel is isomorphic to SLₙ(O): a representative with determinant ε^(nk)
has the unique determinant-one representative ε^(−k)g. Multiplying the lift
by ε^j I changes k by j and leaves this correction unchanged. The map is
multiplicative and SLₙ(O) ∩ ⟨εI⟩ is trivial. The implementation uses all 2n
cosets and rewrites every conjugate of every GL and scalar relator. It uses
componentwise addition, including inverse letters, rather than the XOR
operation which sufficed for C₂ × C₂ in rank two.

In either case the scalar centre of SLₙ(O) is μₙ(O). If ζ generates the
torsion units of order m, this centre has order d=gcd(n,m), generated by
ζ^(m/d)I. That matrix lies in the root finite stabilizer; its verified word
is rewritten into the SL alphabet. PSL is the quotient by this word when
d>1. In odd rank −I is not in SL. In particular, rank-three Eisenstein PSL
requires killing an element of order three, rather than −I. Projective
matrix checks and word certificates use the exact scalar λ with λⁿ=1.

The sections below retain the original rank-two certificate conventions;
their square-unit quotient and {±I} centre are special cases of this argument.

## 1. Arithmetic group and complete vertices

The algebra is M₂(K), its order is M₂(O), and its unit group is GL₂(O).
The lattice is O², with integral basis `(e1,w*e1,e2,w*e2)`.
For D5, w²−w−1=0; for D8, w²−2=0; for D12, w²−3=0;
for D13, w²−w−3=0; and for D40, w²−10=0.
The real space of symmetric forms has
dimension six, and the positive cone consists of forms positive definite
at both real embeddings. Forms have common trace minimum one:

    q_F(v) = Tr_K/Q(vᵗ F v),       F^g = gᵗ F g.

The adapter checks integral unimodularity and commutation with multiplication
by w. Thus it enumerates precisely field-linear lattice automorphisms;
it does not silently enlarge the group to integral semilinear isometries.
Scalar units remain matrices in this group and have not been factored out.

The unchanged M2 validator in `src/Validation.jl` verifies exhaustive minima,
perfection, the complete finite field-linear stabilizers, inequivalence of
representatives, every irredundant facet, and its exact outward neighbour and
transporter. See `BACKEND_AUDIT.md`, [the user guide](USER_GUIDE.md), and the retained M2 results.
The representative counts are two for D5/D8, three for D12, nine for D13
and 52 for D40. Each complete facet
list has a verified neighbour in those representatives. By connectedness of
the perfect-cone adjacency tessellation, this closed neighbour enumeration
contains every perfect-form orbit: starting from a known representative,
follow an adjacency path to any other cone. Induction along that path never
leaves the recorded orbits. Completeness uses this geometry and the verified
local facet lists, rather than a guessed stopping count.

The historical counts 8 (D5) and 7 (D8) count **local directed facet orbits**
under the respective finite form stabilizers. Both fields have four global
unoriented edge orbits. These counts have different meanings.

## 2. Exact positive ridges and complete cycles

`TwoCells.jl` stores evaluation rays in the dual of the form-coefficient basis
C₁,…,C₆. It constructs the nonsingular matrix

    M_ij = Tr_K/Q(tr(C_i C_j)),

and checks `M * coordinates(v*vᵗ) == evaluation_ray(v)` for every supporting
vector. All cone operations use these dual coordinates consistently.

Intersecting every pair of actual facets and retaining rank-four intersections
finds every ridge. Supporting-vector index sets identify faces of this fixed
cone exactly; primitive rational extreme rays are stored as an additional
coordinate certificate. Each ridge belongs to exactly two irredundant facets.
Its supporting vectors have K-rank two, witnessed by a nonzero 2×2 determinant.
The sum B_R = Σ vvᵗ lies in its relative interior and is positive at both
embeddings. Positivity follows because the two independent K-vectors stay
independent at each embedding; exact positive definiteness is also checked.

The absolute ridge stays fixed during a walk. A state stores the actual
normalized form, its representative index and frame a, and both actual
incoming and outgoing facets. With a local stabilizer reduction s and
canonical neighbour F_j^b, crossing changes the frame to **b s a**.
Correspondingly C_(F^g) = g⁻¹ C_F g⁻ᵗ, so the absolute minimal vectors are
obtained with a⁻¹. Every crossing checks the shared facet and ridge in these
coordinates and evaluates its prepended matrix word.

The transverse link of an interior ridge in a locally finite exact fan is a
finite circle. Taking the other incident facet follows that circle once.
The walk accepts only the first return to the complete oriented geometric
state; returning to a representative index is insufficient. An unexpected
repetition or resource limit is an error/partial checkpoint. On closure the
frame h belongs to the initial finite form stabilizer, though h need not be I
and need not fix the initially stored local ridge in its original frame.
The relation is W_walk W_h⁻¹. A finite-stabilizer word for h is recorded and
verified as an actual matrix, including its sign.

Both directions of every ridge of every representative are processed.
Reverse cycles match the same actual forms in reverse order, exchange their
incoming/outgoing facets, and record finite-stabilizer corrections between
their frames. Raw words need not coincide. All changes of starting flag are
covered because every local ridge is processed independently. Given any
global two-cell, send one incident perfect cone to a recorded representative;
its ridge then occurs in that representative's complete local list. Hence
this deliberately redundant list meets every required global two-cell orbit.
We do not claim to have minimized the number of global two-cell orbits.

The coverage ledger has no pending entries. Per-loop caches identify the
field, action, exact mathematical source data and schema. Loads replay all
crossings and compare the entire serialized certificate with its exact
reconstruction. Incomplete, stale and corrupted data cannot certify a result.

## 3. Graph of groups, inversions and coherent lifts

`Presentations.jl` classifies global unoriented edges by actual endpoint
forms. For the edge (F_i,F_j^b), all possible identifications from its first
endpoint to F_i are in G_i, and those from its second endpoint to F_j are
in b⁻¹G_j. Testing these complete finite cosets against both endpoints
exhausts the edge identifications, and retains their matrix witnesses p.

The oriented edge group is E = G_i ∩ b⁻¹G_j b. Its second-end inclusion is
u ↦ b u b⁻¹. If the two endpoints have the same representative, all possible
reversals lie in G_i b; testing this complete coset gives the full setwise
edge stabilizer H. Reversal squares belong to E. In D8 one edge has |E|=2
and |H|=4, and its reversing matrices square to I. There are no inverted
edges in D5. These reversal squares are measured as matrices, not inferred
from the order of H. The implementation also retains order-four elements
with square −I in finite stabilizers; this is checked separately.

Every original edge is subdivided equivariantly. The midpoint vertex group
is its full H, including reversals, rather than E alone. The new half-edges
have no inversions: their original-vertex and midpoint-vertex types differ.
This is an equivariant CW subdivision, with the same underlying contractible
space; subdividing the boundaries of the original two-cells changes no
attaching disk. We need not add higher-dimensional stabilizer data to compute
the fundamental group.

Choose coherent vertex lifts l_v along a spanning tree of the subdivided
quotient graph. Every finite vertex group becomes l_v⁻¹ G_v l_v.
For a half-edge represented by p from original vertex i to midpoint q,
its stable matrix is t_d = l_q⁻¹ p l_i. Tree choices make these matrices
exactly I, so their words can be empty. This is a simultaneous change of
lifts and stabilizers; no nonidentity matrix has been declared trivial.
The inclusions are imposed as

    t_d (l_i⁻¹ u l_i) t_d⁻¹ = l_q⁻¹ p u p⁻¹ l_q.

Finite vertex groups are given genuine presentations. Exact enumeration
provides their elements, multiplication tables, generators and normal-form
words w_x. For every element x and generator g, the relation
`w_x g = w_(xg)` reduces any word to its unique table element. These
relations therefore present the complete finite group, including −I;
the group order is not being used as a substitute for relations.

For completeness of the path rewrite, take an original canonical crossing b,
its source half-edge d with witness p, and the reverse half-edge d′ with
witness p′. Let r in G_j reduce its target facet to the stored reverse facet
orbit. Then z = p′ r b p⁻¹ is in the midpoint stabilizer. In coherent
alphabets put R_j = l_j⁻¹ r l_j and Z = l_q⁻¹ z l_q. The subdivided crossing
word is

    T_b = R_j⁻¹ T_d′⁻¹ Z T_d,       matrix(T_b) = l_j⁻¹ b l_i.

This word describes the two half-edges and their midpoint identification,
not an arbitrary word with the right matrix. The source correction s_src in
G_i is a different matrix: a_next = b s_src a. Put
S_i = l_i⁻¹ s_src l_i. Prepending T_b S_i at each crossing telescopes the
intermediate lifts. On closure we append the inverse
of the initial vertex word for l_i⁻¹ h l_i. Thus the stored ridge relators
are the actual lifted attaching paths required by the group-action theorem.

## 4. Faithfulness, rather than only matrix relations

The positive perfect-cone faces are dual to the cells of the normalized
well-rounded complex. It is contractible; its vertices, edges and two-cells
are dual to the dimension-six cones, dimension-five facets and positive
dimension-four ridges used above. Subdivision preserves this property.

We have all vertex orbits and their full stabilizers, all subdivided edge
orbits with both inclusions, coherent tree lifts, and at least one attaching
loop for every two-cell orbit. The group-action presentation theorem therefore
identifies the raw abstract presentation with GL₂(O). Redundant two-cell
relations follow from the same theorem and do not alter the group. Matrix
evaluation verifies the chosen realization of this isomorphism, rather than
establishing its injectivity on its own. In particular the ineffective central
element −I survives through the finite stabilizer groups.

Equivalently, a graph loop can be filled in the contractible complex; a compact
filling meets only finitely many local pieces. Ridge-cycle relations and
stabilizer identifications account for the resulting moves. Unbounded cones
do not create an additional infinite-word condition. Higher-dimensional
attaching cells do not change this fundamental-group calculation.

## 5. Constructive words and the passage to SL

For a GL matrix g, `matrix_to_word` sets B to the sum of every root cone
generator and Q = g⁻¹ B g⁻ᵗ, which is an interior point of the target cone.
If Q is outside the current cone, cross a separating facet. The exact
objective Tr_K/Q(tr(F_current Q)) strictly decreases. It is positive and
lies in a fixed rational denominator lattice: Q is fixed rational field
data, there are finitely many representative forms, and frames are integral
O-matrices. Thus there cannot be infinitely many decreases. At termination,
interiority and the exact tessellation give the target actual cone uniquely.
The remaining finite stabilizer correction produces a word for g. The
implementation checks each decrease and the final matrix; explicit step/time
limits return errors rather than words claimed to be complete.

Choose ε=w for D5, 1+w for D8/D13, 2+w for D12 and 3+w for D40.
Hecke `unit_group(O; GRH=false)`
verifies O× = {±ε^m}, including that the free exponent of the chosen ε is ±1.
Constructive reduction supplies z with matrix εI. Since the GL presentation
is already proved faithful, adding z=1 gives exactly GL₂(O)/⟨εI⟩.
This is an algebraic central quotient; no trivial scalar action or
contractibility of a quotient complex is assumed.

The determinant modulo squares maps this quotient onto C₂×C₂; its computed
generator images attain all four cosets. Its kernel is isomorphic to SL₂(O):
a lift with determinant ε^(2k) has the unique determinant-one representative
ε^(−k)g in its scalar coset. This rescaling is independent of the lift and
is multiplicative. The intersection SL₂(O) ∩ ⟨εI⟩ is trivial because ε has
infinite order.

`SLExtraction.jl` uses a prefix-closed four-element Schreier transversal.
It rewrites every conjugate of every GL relator and the added scalar relator,
with all four transversal representatives. Each Schreier generator retains
its GL lift word and the scalar correction determined by its exact unit
exponent. The Reidemeister–Schreier theorem therefore gives a presentation
of the full kernel, hence SL₂(O). Every corrected matrix has determinant one,
every relation evaluates to I, and the identified word for −I evaluates
to −I. SL has not been treated as a finite-index subgroup of the original GL.

Finally, `Simplification.jl` uses the OSCAR/GAP Tietze isomorphism, retaining
old-to-new and new-to-old generator words. Its isomorphism status is checked,
both matrix maps commute, and the simplified relators evaluate exactly.
The unsimplified presentations and their relator provenance remain available.
Thus simplification preserves the groups already established above.

## 6. Additional fields, PSL and the class-number-two lattice

The D12 unit has field norm +1. The determinant-square kernel construction
above depends on squares in the complete unit group, not positivity or field
norm, so no change of subgroup is made in that case. PSL is formed only
after the faithful SL presentation, by adding the identified −I word as a
relation. Its saved matrices are SL lifts and its evaluated relators may
equal I or −I. This is a separate central quotient from GL/⟨εI⟩.

The geometry and group-action theorem apply to M₂(O) acting on the one fixed
free lattice O² without a PID assumption. The adapter contains no O-gcd,
principal-ideal normalization, Euclidean division or completion of arbitrary
O-vectors to bases. Its primitive-ray gcd is over rational/integer cone
coordinates and does not impose unimodularity over O. Hecke minima are
exhaustive on the underlying rank-four integer lattice. The Euclidean
row-reduction helper in `Literature.jl` is isolated to the four comparison
fields and is not used for D40 enumeration, words or presentations.

For D40 the ideal 2O+wO has norm two and is nonprincipal: principality would
give x²−10y²=±2, impossible modulo five. The vector (2,w), with integer
coordinates (2,0,0,1), is integer-primitive but not O-unimodular. The bounded
enumeration regression includes it at trace value 28 for F=I. This test is
about a bounded enumeration, not a claim that it is a minimum of I. The
ordinary class group is computed with GRH=false and has order two. No
uncomputed narrow class group is substituted for it, and no other rank-two
projective lattice is required to compute GL₂(O) of the stated free lattice.
