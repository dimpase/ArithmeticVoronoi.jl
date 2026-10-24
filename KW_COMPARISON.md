# Kirchheimer–Wolfart comparison

The comparisons concern the PSL theorems (9.1), (10.1), (11.3) and (12.2)
of F. Kirchheimer and J. Wolfart, *Explizite Präsentation gewisser Hilbertscher
Modulgruppen durch Erzeugende und Relationen*, J. reine angew. Math. 315
(1980), 139–173, and the SL central extension in section 13. D12 uses (11.3),
not the index-two projective overgroup in (11.2).

The local sources read were the user's `KWeng.tex` and `KWger.tex` in
`/home/dima/work/papers/sl2translation/`; identical copies are in `papers/sl2`.
Source hashes are retained in each comparison record. The local PDFs are
compiled transcriptions, not an independently inspected original scan.

The comparison is outside the computation package. Load it explicitly:

```julia
include("validation/KirchheimerWolfart.jl")
using .KirchheimerWolfart
compare_kw(R)
```

For a fresh checkout, `compare_kw(R)` uses the previously recorded hashes in
`references/kw_sources.json` and labels them `recorded_hashes_only`. It does
not require the original account's directories. To check your own copies,
use `compare_kw(R; reference_source_directory="/path/to/transcriptions")`;
both files must match the recorded hashes, and the result is labelled
`checked_local_sources`. The TeX transcriptions are not redistributed here.

## Reference alphabets and signs

`validation/KirchheimerWolfart.jl` constructs the four PSL reference presentations in the common
alphabet U=τ(1), V=τ(e), H=diag(e,e⁻¹), J=[0,1;−1,0]. It generates both
conjugation relations directly from exact multiplication by e² in the basis
(1,e). The coefficient pairs (α,β), for e²=α+βe, are (1,1), (1,2), (−1,4)
and (1,3). All torsion relations are evaluated with the displayed J, which
is the negative of the other common choice [0,−1;1,0].

For SL the fifth generator Z maps to −I. It is central of order two,
J²=(HJ)²=Z, while (JU)³=(HJV)³=I. The extra D12 relation has square Z.
Quotienting by Z gives the published PSL presentation. Hence the kernel of
the SL reference's surjection onto SL lies in the central subgroup generated
by Z; Z has nonidentity matrix −I, so the kernel is trivial. The explicit
reverse substitutions also establish that the reference matrices generate
the computed SL group. This proves the finite SL reference is faithful,
conditional on the published PSL theorem.

## D8 transcription discrepancy

English (10.1) ends the third relator with τ(e)⁻². The German transcription
has τ(e)⁻⁵, and the introduction to English section 10 gives k=5. Since
e=1+√2 satisfies e³=2+5e, the correct relation is
H V H⁻¹ U⁻² V⁻⁵=I. The English displayed word instead equals τ(3e), which
is neither I nor −I. The original English translation has not been edited.
This is documented as a discrepancy in the supplied transcription; no error
in the published paper is asserted without inspecting its scan.

## Meaning of the comparison certificates

`results/fields/D*/KW/comparison.json` records substitutions in both
directions for D5, D8, D12 and D13. `reference.json` contains the exact
reference matrices and SL/PSL words. Each reference-to-computed SL generator
has a replayable Voronoi path certificate. The opposite direction uses exact
row operations with strictly decreasing absolute norms of lower entries;
each actual quotient is retained. This routine is restricted to the four
comparison fields and is never invoked by the generic Voronoi engine or D40.
Its explicit search-radius/step bounds can fail, rather than assuming an
unproved global Euclidean algorithm.

There is no field conjugation or matrix conjugation: our basis is (1,w),
and e has the exact coordinates recorded in `field.json`. Both round trips
preserve the actual SL matrices, including their central signs. Substituted
relators and abelianizations are checked.

The proof route is **independent faithfulness of both matrix presentations**:
the computed group is faithful by geometric coverage and the group-action
theorem; the reference is faithful by the published PSL theorem and the
central-kernel argument above. Consequently the commuting matrix
substitutions define mutually inverse abstract group maps. We do not claim
to save formal relator derivations or a Tietze chain connecting the reference
alphabet directly to the computed alphabet. Tietze maps within each computed
presentation are retained separately.

The D8 hand-checkable fixture from the handoff is tested in the preserved
original simplified generators a,b, including b⁴=−I, the translations τ(1)
and τ(w), the diagonal H, and the sign of J. The D5/D8 original outputs and
their maps remain in their original directories and in the baseline archive.
D40 has no Kirchheimer–Wolfart comparison claim.

The PSL substitutions use only the four reference generators U,V,H,J.
The auxiliary Z letters from an exact SL row reduction are deleted upon
projection, and their discarded central signs are recorded separately.
Exact SL-lift words in the five-generator alphabet are retained as well.
Both projective round trips and substituted relators are checked modulo ±I.
