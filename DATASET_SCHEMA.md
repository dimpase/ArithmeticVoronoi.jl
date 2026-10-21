# Arithmetic-group dataset, schema v1

The local export is `results/dataset/v1/dataset.json`. It describes arithmetic
objects independently of a particular presentation. It is a proposal format,
not an existing or accepted LMFDB table. LMFDB labels and generator maps are
null until verified; the attempted D5/D8 field pages returned browser challenges.

## Encoding and identities

The root `schema` is `arithmetic_groups_dataset_v1`. Each field records its
discriminant, defining polynomial in ascending coefficients, maximal order,
integral basis (1,w), and coordinate basis (e1,w*e1,e2,w*e2). Field elements
are arrays of **rational strings** in ascending powers of w. Matrices are
arrays of rows of those elements. No decimal embedding or display string is
used for mathematical arithmetic. Counts and word indices are JSON integers;
timing data is explicitly nonmathematical floating metadata.

A word is an array of signed, one-based generator indices; the empty word
is I and products are evaluated in written order. A negative index is the
inverse matrix. Action convention: F^g=gᵗFg on columns, and frames update
by b*s_src*a. All matrices are K-linear on the free lattice O².

An `object_id`, e.g. `real_quadratic:D12:Omax:free_rank2:K_linear:SL`,
identifies the field/order/lattice/group specification. The `arithmetic`
object gives determinant and scalar-quotient conditions explicitly.
Presentation versions are SHA-256 identifiers of raw presentation artifacts;
computation-run identifiers separately hash the recorded run. Different valid
presentations may have the same arithmetic object identifier.

GL has unit determinant and no scalar quotient, and includes all scalar
units. SL has determinant one and retains −I. PSL is SL/{I,−I}; its matrices
are determinant-one lifts, and its relators may evaluate to either I or −I.
Raw PSL relators retain their measured lift signs. The intermediate
GL/⟨eI⟩ in `group_maps.json` retains −I and is explicitly **not PGL**.
Its determinant-square quotient C₂×C₂ uses unit-group square classes,
independently of total positivity or the sign of the field norm of e.

## Artifacts and maps

Every artifact reference has `path`, `sha256` and `bytes`, relative to the
bundle root. The root artifact index includes every supplied file. Hashes
identify data; they do not prove its mathematical correctness.

Each completed field has exact field/unit/ordinary-class-group metadata,
independent stage statuses, run timings, raw GL/SL/PSL presentations and
simplifications, and abelianizations. The simplification files retain both
generator maps. `group_maps.json` records the fundamental scalar word,
the −I words, the central scalar quotient, determinant-square images and
Schreier data, corrected SL generator lifts, the SL-to-GL embedding words,
and the SL-to-PSL quotient. Domains and codomains are named explicitly.

`checkpoint.json` contains the versioned exact exporter payload and its hash.
`two_cells/coverage.json` and `two_cells/loops/` contain all local ridges in
both orientations and exact reverse-frame corrections. Geometry counts
distinguish perfect-form orbits, local directed facets, global unoriented
edges, local/local-stabilizer ridge counts, and oriented attaching loops.
Global two-cell orbit minimization is not performed and its count is null.
Constructive-word samples retain paths, frames, crossing labels, objective
values, a common objective denominator, the root-stabilizer correction,
alphabet fingerprint, and the central lift sign. They are not normal forms.

Kirchheimer–Wolfart records include both substitutions, field/sign conventions,
source hashes, the D8 discrepancy, and the independent-faithfulness proof
method. D40's literature-comparison status is not applicable, not verified.

The pinned Julia environment, exact package sources and checker are copied
under `software/`. Provenance gives dependency versions, the pinned upstream
backend commit, code hashes and input/output hashes. The workspace package
is not attributed to a dedicated tracked commit: `workspace_git_commit` is
null and the exact source artifact hashes identify the implementation.

## Verification entry points and strength

From the repository root, after `git lfs pull` and unpacking the dataset as
shown in [the user guide](USER_GUIDE.md#optional-historical-archives), install this checkout's Julia dependencies and run:

```sh
julia --startup-file=no --project=. \
  -e 'using Pkg; Pkg.instantiate(;allow_autoprecomp=false)'
julia --startup-file=no --project=. scripts/verify_dataset.jl \
  results/dataset/arithmetic_groups_dataset_v1

julia --startup-file=no --project=. scripts/verify_dataset.jl \
  results/dataset/arithmetic_groups_dataset_v1 --complete
```

To use the software pinned inside the bundle instead of this checkout:

```sh
julia --startup-file=no --project=results/dataset/arithmetic_groups_dataset_v1/software \
  -e 'using Pkg; Pkg.instantiate(;allow_autoprecomp=false)'
julia --startup-file=no --project=results/dataset/arithmetic_groups_dataset_v1/software \
  results/dataset/arithmetic_groups_dataset_v1/software/scripts/verify_dataset.jl \
  results/dataset/arithmetic_groups_dataset_v1 --complete
```

Dependencies must be installed from the included pinned Manifest; no backend
executable is needed to validate the exported exact enumeration payload.
Nothing in the interchange format is evaluated as Julia, GAP or shell code.

The default **algebraic** level reconstructs exact fields and matrices,
checks membership, raw and simplified relators with the right central
semantics, both matrix maps, unit/class data, abelianizations and quotient
corrections, then replays all supplied word certificates. It also reuses
the existing enumeration validator to construct the fan context. It does
not re-certify attaching-loop coverage or the abstract saved Tietze maps.

The **complete** level additionally replays every oriented ridge cycle and
coverage witness, reconstructs the group-action presentation and the SL/PSL
extractions, and repeats the Tietze isomorphisms. A mismatch in the saved
Tietze maps is reported as a failure to re-certify, never ignored on the
basis of matrix equality. The interpretation as a faithful arithmetic-group
presentation invokes `PRESENTATION_JUSTIFICATION.md` and its theorem
hypotheses. Neither level is a proof-assistant formalization.
