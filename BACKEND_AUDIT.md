# Arithmetic Voronoi backend acceptance gate

Executed 2026-10-17 for `oscar_arithmetic_voronoi_codex_handoff.md`.

This is a historical audit. The C++ sources and build now live in
[PerfectForms.jl](https://github.com/dimpase/PerfectForms.jl); paths below
describe the original workspace. Use [the user guide](USER_GUIDE.md) for
current installation and native Julia computations.

**Verdict: usable after the listed local fixes.** This verdict covers the
GMP serial perfect-form/facet-neighbour path for the free rank-two lattices
over Z and the maximal orders of discriminants 5 and 8. It does not certify
upstream's full complex, presentation routines, other arithmetic types, or
larger dimensions. No broad arithmetic-engine repair was needed. Keep the
replaceable process boundary and independent validator.

## What actually runs

The isolated Julia package implements M1 and M2, using a small separately
maintained C++ exporter around the pinned upstream enumeration callbacks.
The exporter captures fresh neighbours, stabilizers and transporters before
they are lost by a top-level query; it constructs no lower complex.

| Example | Perfect orbits | Minimal-vector pairs per representative | Stabilizer orders | Actual facets | Labelled facet orbits |
| --- | ---: | --- | --- | --- | ---: |
| GL(2,Z) | 1 | 3 | 12 | 3 | 1 |
| GL(2,Z[w]), w²−w−1=0 | 2 | 10, 6 | 20, 12 | 22, 6 | 8 (6+2) |
| GL(2,Z[w]), w²−2=0 | 2 | 6, 12 | 12, 16 | 6, 24 | 7 (2+5) |

All required facets close, and all representatives have exact minimum one.
Counts are measured output, not hardcoded enumeration bounds. The final
machine-readable records are [D5](results/D5.json), [D8](results/D8.json),
[classical2](results/classical2.json) and [summary](results/examples_summary.json).

## Pinned build and ownership

- Top source: `MathieuDutSik/polyhedral_common` commit
  `cfcd79dbddf67a494a1c71b776174b190da5209a`.
- All seven recursive submodule pins and the nauty fork pin are in
  the former `backend/pins.json`. Submodules were initialized without
  `--remote`. The selected nauty commit is
  `cb9dc822f9d1aa3a02b1dcb7a18dd9609d838b1e`.
- Source and build are separate directories under
  `/tmp/hilbert-ford/arithmetic_voronoi_backend/`. Only the top repository's
  `CMakeLists.txt` is patched; the symlinked source directories and submodule
  ownership were inspected. No upstream arithmetic headers were edited.
- Release: GCC 16.2.0, C++20, CMake 4.3.4/Ninja, `-O1 -DNDEBUG`,
  `-DSANITY_CHECK`; one build job, 10 GiB virtual-memory ceiling. GMP 6.3.0,
  Eigen 3.4.1, Boost 1.90. MPI, FLINT, SCIP, Boost number types and GMP pool
  disabled. GMP integers/rationals are the only exporter arithmetic.
- Julia 1.13.2, OSCAR 1.8.2, Hecke 0.39.22, GAP.jl 0.17.6 / GAP 4.16.1,
  JSON3 1.14.3. The local [Project](Project.toml) and [Manifest](Manifest.toml)
  pin the environment. Resolving this isolated environment adjusted bundled
  networking JLLs for installed Julia; the root projects were not edited.
- [Build records](results/build/) retain the command, compiler/dependency
  versions, CMake cache, source diff, recursive pins, compile commands, nauty
  configure log, warnings/logs, source SHA and executable SHA. Historical
  refreshes are separate from fresh compilation timings.

Nine pristine source trees, including recursive submodules and nauty with
their license notices, are retained in [pinned_sources.tar.gz](backend/pinned_sources.tar.gz).
Its internal manifest hashes each nested source tar; the outer SHA is in
[pinned_sources.sha256](backend/pinned_sources.sha256). The 8.5 MB archive was
created and independently read back with `scripts/archive_backend_sources.py`.
Local changes remain separately reviewable; this is source preservation,
not an offline dependency installer.

The accepted release SHA-256 is
`36af80176b246cd47a08e295158c8842c292cf10fc211333d161aeb541a8400c`.
Its explicit [configuration](backend/config.json) selects the retained
`backend/bin/AV_ExportPerfectOrbits`, rather than any newest build or PATH entry.
The original driver SHA is
`b9cea76e9a2f54f092a38e0d68869fc5a622069b51442f59904bd780aba23c2b`.
The sanitizer SHA is
`ac10134c08fa7bed7debdd14b6baed081a4b6411d0d0fb269d5b5452e6fb7b74`.

The official target's first serial build took 290.49 seconds including nauty.
The sanitizer compilation took 163.02 seconds. The latest release build record
is an incremental no-op/configuration refresh; its subsecond duration must not
be reported as a fresh compile. The sanitizer `-Wall -Wextra` build emitted no
diagnostics. The retained official-target/nauty build log also has no compiler
warning diagnostics. This records the selected builds, not a generic guarantee
about every target or compiler.

## Findings and local fixes

| Pinned source observation | Executed evidence / accepted treatment |
| --- | --- |
| `ComputeComplex` is declared, but the serial driver calls full-complex construction regardless. | With `ComputeComplex=F`, the official classical probe exported counts `[1,1,0]`; D5 exported `[2,4,5,4,1,0]`. These include lower levels. Inputs, output, status and timings remain in [upstream_probe](results/upstream_probe/). Use the narrow exporter, not the flag as a performance guarantee. |
| GAP `__PERFCOMP_Write_t_space` has cache emission commented out. | The adapter writes its exact namelist and owns all cache configuration. Backend archives are disabled with `CacheFile="none"`. |
| GAP `PERFCOMP_general_query` concatenates shell text and reads GAP code, without explicit child-exit validation there. | Julia uses an argv vector, private files, captured stdout/stderr, process status/signal and a request nonce. JSON contains exact rational strings; no GAP evaluation, `eval`, or `ReadAsFunction` is used. The executable and request paths with spaces pass. |
| `GetBinaryFilename` can choose the newest matching build. | One absolute real path and expected SHA are configured. Replacing the executable is rejected; sanitizer builds stay separate. |
| Upstream Boost cache output is direct, with limited compatibility metadata. | The adapter owns a versioned field/order/basis/group/options/build/executable key, raw-payload checksum and atomic final rename after mathematical validation. Warm cache loads repeat validation. Truncation and incompatible envelopes are rejected. |
| Full complex prints the result of `are_product_zeros`, then returns even if it is false. | No boundary/complex data is accepted. Boundary requests fail explicitly; a cache carrying `boundary_square=false` is rejected despite valid JSON/checksum. Actual bad boundary products were not observed or synthesized inside upstream: this is an avoided control-flow risk, not a claim to have repaired the full-complex path. |
| Real-quadratic test generation exists, but the inspected invocation is commented out; workflow is scheduled. | The actual D5/D8 constrained computations and independent checks were run locally. No inference of coverage from CI status. |
| `update.sh` uses `--remote`; CMake nauty used a branch name. | Never invoke that updater. The small CMake patch pins nauty and adds the external exporter source target. Configure checks the top commit, recursive submodules, and exact tracked diff. |
| Serial enumeration can return a partial database at an internal time limit. | The exporter disables that internal cap, checks each orbit was processed, and writes completion only after queue exhaustion. A Julia outer watchdog terminates overruns; no success cache is published. |
| The code's RNG is its own `mt19937_64` engine. | `set_random_seed(1)` is called, as well as the C seeds. Merely calling `srand(1)` was insufficient. Final accepted runs and sanitizer use this corrected exporter. |

Relevant pinned source is
[serial driver](https://github.com/MathieuDutSik/polyhedral_common/blob/cfcd79dbddf67a494a1c71b776174b190da5209a/src_perfect/PERF_SerialPerfectComputation.cpp),
[perfect_tspace](https://github.com/MathieuDutSik/polyhedral_common/blob/cfcd79dbddf67a494a1c71b776174b190da5209a/src_perfect/perfect_tspace.h),
[complex](https://github.com/MathieuDutSik/polyhedral_common/blob/cfcd79dbddf67a494a1c71b776174b190da5209a/src_perfect/perfect_complex.h),
and [GAP access](https://github.com/MathieuDutSik/polyhedral_common/blob/cfcd79dbddf67a494a1c71b776174b190da5209a/CI_tests/access_points.g).
The retained patch (`patches/polyhedral_common.patch` in the original layout)
is CMake-only; the exporter (`backend/export_perfect.cpp`) is approximately
220 lines. The useful
upstream strengths are exact GMP arithmetic, constrained field-linear tests,
existing neighbour/transport records and finite-group computations. The
successful narrow gate is not a whole-repository engineering endorsement.

## Independent mathematical checks

[Model tests](test/model.jl) check both maximal-order bases, exact trace
evaluation, dimension six, positivity, encode/decode round trips,
`A_F E = Eᵗ A_F`, action composition and matrix conversion. A genuine integral
Galois isometry of the identity trace form is rejected for not commuting with
field multiplication. D8's nonsymmetric multiplication matrix catches the
otherwise easy transpose error. Nonidentity lattice bases and rational Gram
scaling test Hecke's lattice-coordinate and squared-length conventions.

[Validation](src/Validation.jl) independently recomputes exhaustive signed
shortest vectors in Hecke, relative perfection, and exact cone facets in
OSCAR/polymake. Every stabilizer generator is integral, unimodular and K-linear,
and preserves its form. Generated closure equals the independently computed
Hecke number-field lattice automorphism order, proving completeness; `-I`
is retained. Number-field isometry uses the multiple restricted-scalar forms,
not an arbitrary single-Gram GL4(Z) witness.

Representatives are independently inequivalent. Every facet orbit has a
distinct perfect endpoint in its exact outward direction, correct common
minimum face, and transporter satisfying `H=gᵗF_target g`. Full minimum sets
and signed facet identifications map correctly. All exact facets are covered
under the complete stabilizers, with consistent reduction words/permutations;
reverse orbit incidences close. No unresolved equivalence is treated as false.

The global inference is the arithmetic perfect-form tessellation and neighbour
connectedness, specializing the algebra/lattice to `M2(K)` and `O^2`. Exhausting
all certified facets from a perfect seed reaches all orbits; minimum one fixes
the common scale. See Lemma 2.3 and Theorem 2.4 of
[Coulangeon--Nebe--Braun--Schoennenbeck](https://arxiv.org/html/1407.6234).
This does not use the four-dimensional determinant-one H×H slice, quotient
minimal vectors by infinite units, or enlarge the group to semilinear maps.

## Engineering regression results

[Final test transcript](results/test_summary.txt): **231 passing assertions**
(97 model, 6 shortest-vector conventions, 128 process/cache/label contracts).
[Structured contract evidence](results/contracts_summary.json) identifies
the preserved input/log/process directory for each failure case.

- Missing/nonexecutable and relative executable; changed executable SHA.
- Real malformed namelist, unsupported arithmetic and boundary request;
  unwritable output; nonzero child exit; signal termination; timeout;
  missing, truncated, stale-token and incomplete output.
- Executable and request paths containing spaces; repeat computation with
  seed 1 and independent K-linear orbit matching.
- Wrong field, subgroup, order basis, commuting options, schema and build;
  incomplete/truncated/checksum-corrupt cache and missing envelope field.
- Checksum-consistent mathematical tampering: incorrect minimum, subgroup-only
  stabilizer, unresolved facet, invalid transporter and false boundary assertion.
- 81-digit integral matrix coefficients survive exact wire encode/decode.
- An actual backend deadline kills a run. No completed cache appears. Restart
  uses the same problem/build key, recomputes successfully and then warm-loads.

Every completed backend orbit appends a partial progress journal. This initial
adapter does **not** deserialize a partial upstream queue: restart recomputes
from scratch, preserving diagnostics and all unresolved work by refusing to
promote the journal. Only final validated checkpoints are reusable. This is
an explicit limitation for long future enumerations.

Malformed/truncated fresh output is never parsed as executable code, and a
normal process exit alone cannot mark enumeration complete. The process
report distinguishes exported-but-unvalidated from validated-complete output.
Failed records remain diagnostic; none is an accepted mathematical checkpoint.

## Sanitizers, conversions, and timings

The Clang 22.1.8 build uses ASan+UBSan, C assertions (`-UNDEBUG`),
`SANITY_CHECK`, `-Wall -Wextra`, and frame pointers. Classical2, D5 and D8
complete cleanly. Malformed input and unsupported arithmetic/boundary requests
also terminate with their intended statuses and no sanitizer diagnostics.
[Sanitizer report](results/sanitizer/summary.json) includes each input hash,
argv, exit, timing, sampled maximum resident memory, and options. Their outputs
were independently revalidated and matched to release orbits in OSCAR/Hecke.

LeakSanitizer cannot operate in the tracing sandbox; the accepted checks ran
natively. Runtime uses a sampled 10 GiB RSS ceiling and 60 second deadline,
because ASan needs huge virtual shadow mappings. This is not a hard RSS cgroup
limit. The pinned nauty static library and linked GMP/Boost libraries were not
instrumented, so coverage is limited to the exporter/instantiated C++ path.

Selected arithmetic remains `mpz/mpq`; all boundary arithmetic is encoded as
strings and decoded to arbitrary precision. Backend permutation indices use
`uint32_t` and some matrix/orbit dimensions use `int`: accepted examples have
at most 24 signed minimal vectors per cone and two perfect orbits. Larger
cases require revisiting those limits; this audit does not prove generic
overflow safety. The release has ordinary `NDEBUG`; the sanitizer build removes
it. No disabled assertion is being counted as executed release validation.

[Example timings](results/examples_summary.json) separate process and exact
validation time, and cold versus warm caches. First in-process use includes
Julia compilation. The final D5 run used about 0.063 s in C++,
7.64 s in first D5 validation and 0.034 s on a subsequent validated warm load.
These
are small-case timings on this host, not scaling claims. Enumeration and
builds are serial; ten-thread HSR settings are unrelated to this backend.

## Reproduction and next gate

Current commands and prerequisites are in [the user guide](USER_GUIDE.md#optional-c-backend).
The commands below record the historical workspace layout. From the
workspace root, the principal checks actually executed are:

```sh
python3 experiments/arithmetic_voronoi/scripts/configure_backend.py
python3 experiments/arithmetic_voronoi/scripts/build_backend.py --target AV_ExportPerfectOrbits --configure-default
julia --startup-file=no --project=experiments/arithmetic_voronoi experiments/arithmetic_voronoi/scripts/run_examples.jl
julia --startup-file=no --project=experiments/arithmetic_voronoi experiments/arithmetic_voronoi/test/runtests.jl
python3 experiments/arithmetic_voronoi/scripts/build_backend.py --target AV_ExportPerfectOrbits --sanitizer
python3 experiments/arithmetic_voronoi/scripts/run_sanitizer.py
julia --startup-file=no --project=experiments/arithmetic_voronoi experiments/arithmetic_voronoi/scripts/validate_sanitizer.jl
```

M3 remains: recover ordered cycles around codimension-two perfect-cone faces,
with edge stabilizers, inclusions, inversions and transporters. Establish the
equivariant two-skeleton and simple-connectedness before applying Brown's
presentation theorem. A chain boundary or adjacency graph alone is insufficient.
The full-complex boundary failure propagation must be fixed before that path
could be accepted. Constructive matrix-to-word reduction and the scalar-unit
word are also pending. M4 must handle the infinite determinant unit quotient,
then the central quotient and index-four kernel with tracked SL matrices.
**No GL/SL presentation or global word solver is implemented or claimed.**

Upstream attribution and license texts are retained in
the original `backend/NOTICE.txt`. No Magma runtime or borrowed
SLnCohomology code was used. The acceptance gate was completed without
publication, upstream pushes or maintainer contact. Subsequently, the user
requested three upstream issue reports; exact texts, reproductions and verified
receipts are in [upstream_reports/README.txt](upstream_reports/README.txt).
