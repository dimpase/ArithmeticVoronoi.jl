Arithmetic Voronoi computations through OSCAR
============================================

Start with USER_GUIDE.txt for installation, reading the supplied data,
running a field, and converting a matrix to a verified generator word.
All commands in that guide run from this package directory. The sections
below are the original engineering and audit notes; their machine-specific
paths describe the accepted experiment, not a portable installation.

Implemented and executed on 2026-10-17. This isolated package completes the
M0 engineering gate, M1 exact model, M2 labelled perfect-form enumeration,
and the two-cell, GL presentation and SL extraction follow-up. The original
backend and arithmetic validator are preserved. BACKEND_AUDIT.md records the
backend acceptance decision; PRESENTATION_JUSTIFICATION.md explains why the
new presentations are complete. Existing Ford/HSR experiments have their own
projects and are not dependencies of this package.

Exact completed examples
------------------------

Group / order                  Perfect orbits  Stabilizer orders  Facet orbits
GL(2,Z)                              1              12                  1
GL(2,Z[w]), w^2-w-1=0, D=5           2              20,12               8
GL(2,Z[w]), w^2-2=0,   D=8           2              12,16               7

For D=5 the two cones have 10 and 6 minimal-vector pairs, 22 and 6 actual
facets, and 6 and 2 facet orbits. For D=8 these numbers are respectively
6 and 12, 6 and 24, and 2 and 5. All facets are accounted for; pending_facets
is empty. The orbit counts are discovered, then validated, not stopping rules.
The original 231 regression assertions and all six ASan/UBSan runs passed.
Sanitizer output also passed the independent OSCAR/Hecke mathematical validator.
The extended regression results are in results/two_cell_tests.txt and
results/two_cell_test_summary.json; the backend was not changed for this stage.

The facet-orbit column counts local directed facet orbits, not global
unoriented edges. Each quadratic example has four global unoriented edges.

Run on this machine
-------------------

From /home/dima/work/hilbert-ford:

  julia --startup-file=no --project=experiments/arithmetic_voronoi \
    experiments/arithmetic_voronoi/scripts/run_examples.jl

  julia --startup-file=no --project=experiments/arithmetic_voronoi \
    experiments/arithmetic_voronoi/scripts/run_two_cells.jl

  julia --startup-file=no --project=experiments/arithmetic_voronoi \
    experiments/arithmetic_voronoi/test/runtests.jl

run_two_cells.jl accepts optional discriminants, e.g. 5 or 5 8, and defaults
to both cases. It reloads and validates the M2 caches, resumes completed
per-ridge loops, rebuilds GL/SL presentations, and simplifies with tracked
Tietze maps. All computations use exact arithmetic. Ordinary execution uses
one Julia computation thread; no ten-thread speedup is claimed here.

The accepted release executable is retained at backend/bin/AV_ExportPerfectOrbits,
so normal runs do not depend on /tmp. backend/config.json selects its absolute
path and SHA-256. Changing that binary invalidates its configuration and caches.
It is a local Linux build linked to this machine's libraries; rebuild on a
different installation. No executable is selected from PATH or by build date.

Rebuild from pinned source
--------------------------

Prerequisites: Julia 1.13, a C++20 compiler, CMake, Ninja, Git, Python 3, GMP,
Eigen, Boost.Serialization, and Linux prlimit. The accepted release uses GCC
16.2.0, CMake 4.3.4, GMP 6.3.0, Eigen 3.4.1 and Boost 1.90. The sanitizer build
uses Clang 22.1.8. Julia packages are pinned by Project.toml and Manifest.toml:
OSCAR 1.8.2, Hecke 0.39.22, GAP.jl 0.17.6 (GAP 4.16.1), JSON3 1.14.3.

These commands have been executed with this project:

  julia --startup-file=no --project=experiments/arithmetic_voronoi \
    -e 'using Pkg; Pkg.instantiate(;allow_autoprecomp=false)'
  python3 experiments/arithmetic_voronoi/scripts/configure_backend.py
  python3 experiments/arithmetic_voronoi/scripts/build_backend.py \
    --target AV_ExportPerfectOrbits --configure-default

configure_backend.py clones when necessary, refuses a different top commit,
initializes exact recorded submodule commits without --remote, applies the
small CMake patch, and checks the tracked source diff. Defaults are:
  /tmp/hilbert-ford/arithmetic_voronoi_backend/polyhedral_common
  /tmp/hilbert-ford/arithmetic_voronoi_backend/build
Use matching --source and --build options on both scripts to relocate these.
The nauty fork is pinned by the CMake patch. The build uses one job, -O1,
SANITY_CHECK, and a 10 GiB virtual-memory ceiling. Optional MPI, FLINT, SCIP,
Boost number types and the GMP pool are disabled. Only GMP mpq/mpz arithmetic
is accepted by the exporter. Build logs, flags, dependency versions, patches,
compile commands and executable hashes are retained in results/build/.
--configure-default explicitly installs the release binary into this project;
a sanitizer build cannot become the default by accident.

The upstream source trees are external dependencies, fetched at their recorded
commits by configure_backend.py. This repository distributes our exporter,
the small CMake patch, dependency pins and reference license texts; it does
not vendor an upstream source snapshot or distribute a backend executable.
BACKEND_AUDIT.md describes the original audit, which also made a local source-
preservation archive. That archive and its helper are now ignored local files
and are absent from the published Git tree and its LFS archive set.

To reproduce the separate official-driver audit, build:

  python3 experiments/arithmetic_voronoi/scripts/build_backend.py \
    --target PERF_SerialPerfectComputation

The exact argv, input.nml, stdout/stderr, exit status and timing for both
executed official probes are in results/upstream_probe/{classical2,D5}/.
They demonstrate that ComputeComplex=F still constructs lower levels. The
adapter uses the separate exporter instead of this full-complex driver.
Those .g files are diagnostic text; the adapter never evaluates them.

Sanity/sanitizer reproduction
-----------------------------

  python3 experiments/arithmetic_voronoi/scripts/build_backend.py \
    --target AV_ExportPerfectOrbits --sanitizer
  python3 experiments/arithmetic_voronoi/scripts/run_sanitizer.py
  julia --startup-file=no --project=experiments/arithmetic_voronoi \
    experiments/arithmetic_voronoi/scripts/validate_sanitizer.jl

LeakSanitizer needs native execution outside a tracing sandbox. Runtime uses
a sampled 10 GiB resident-memory watchdog and 60 second deadline: ASan's
shadow mappings require a large virtual address space, so RLIMIT_AS is used
for compilation and ordinary execution, not for ASan runtime. The exporter
and instantiated upstream C++ headers are instrumented; linked GMP, Boost
and nauty libraries are not. Results and coverage are in results/sanitizer/.

API and conventions
--------------------

  using ArithmeticVoronoi
  P = arithmetic_problem(5)       # field discriminant, not radicand
  E = enumerate_perfect_forms(P)  # revalidates compatible warm cache
  E.completion_status            # :complete only after all exact checks
  E.representatives              # PerfectFormRecord objects
  E.neighbours                   # NeighbourRecord objects
  write_enumeration("D5.json", E)

classical_problem() is the GL(2,Z) convention check. arithmetic_problem(D)
supports discriminants 5, 8, 12, 13 and 40, always on the free maximal-order
lattice O^2. D12 uses w^2=3, D13 uses w^2-w-3=0, and D40 uses w^2=10.
backend=:polyhedral_common is implemented; the mathematical
records are independent of that backend. A native backend is a future option.
No Magma installation is required.

The quadratic integral basis is (1,w); column coordinates are
(e1,w*e1,e2,w*e2). Forms lie in the full six-dimensional trace space, with
one common minimum normalization, q_F(v)=Tr(v^t F v) and F^g=g^t F g.
Backend grouped row coordinates are converted with permutation [1,3,2,4]
and transposition of action matrices. D=8 tests detect accidental transposes.
Every group matrix commutes with multiplication by w, is integral unimodular
over Z, decodes to a matrix over O, and has unit determinant over K. Galois
semilinear Z-isometries are explicitly rejected. Hecke uses the number-field
lattice isometry/automorphism API with its restricted-scalar forms.

Each representative stores minimum 1, a complete set of minimal vectors up
to +/- (first nonzero coordinate positive), exact evaluation rays, stabilizer
generators/order, and every exact cone facet. -I is retained in stabilizers.
The rays live in the dual of the specified form-coefficient basis; these are
relative perfect cones, not unconstrained rank-four perfect-form cones.

Each neighbour H and target F_j carry a matrix g with H=g^t F_j g. Thus g
maps neighbour minimal vectors to target minimal vectors. signed_target_facet_map
records signs and target indices. Facet reduction matrices, stabilizer words
and induced permutations all map that facet to its recorded orbit representative.
Word products are in written order; a negative index means inverse. These
are words in the stored finite stabilizer generators, not global GL words.

Validation, data, and restart policy
------------------------------------

The independent validator checks total positivity, exhaustive Hecke minima,
evaluation rank, complete field-linear stabilizers (generated closure and
independent order), distinct arithmetic orbits, every OSCAR cone facet,
outward neighbour directions, endpoint minima, transporter congruence/full
minimum-set mapping, facet-orbit coverage and reverse orbit identifications.
The completion argument uses the connected locally finite arithmetic perfect-
form neighbour tessellation, as stated in Theorem 2.4 and Lemma 2.3 of
Coulangeon--Nebe--Braun--Schoennenbeck, arXiv:1407.6234. This is a theorem-based
closure certificate, not a claim that a finite search bound proves completeness.

results/{classical2,D5,D8}.json contain exact field coefficients, trace Grams,
minimal vectors, stabilizers, cones/facets, reduction words, neighbours and
completion provenance. All arithmetic entries are integer/rational strings;
field entries are coefficient lists in the fixed basis, never decimal embeddings.
results/examples_summary.json separates backend, validation, cold and warm
in-process timings. First use includes compilation; these are not Julia
startup or filesystem-cold benchmarks. results/test_summary.txt and
results/contracts_summary.json record the regression results.

Requests have private directories (deliberately with spaces), exact input,
stdout/stderr, nonce, hashes and process.json. Missing/changed executables,
invalid input, unsupported arithmetic, output failure, child exit/signal,
timeout, missing/truncated/stale output, incompatible caches and invalid
arithmetic have distinct errors. No shell command concatenation or eval is
used for interchange. Only independently validated complete output is
atomically published as a cache. Warm loads repeat the mathematical checks.

Each completed backend orbit appends a PARTIAL progress journal. It is not
a reusable enumeration database. An interrupted computation retains its
input/logs/journal and restarts enumeration from scratch with the same key;
only a fully validated final cache is reusable. This initial implementation
does not resume a half-processed queue. The killed-run restart is tested.
Historical failed tests/probes are retained; accepted current examples and
the final summaries identify the successful artifacts. export_probe/ predates
the final RNG seeding fix and is not an accepted completion certificate.

Completed two-cells and presentations
-------------------------------------

All groups below are K-linear on O^2. GL includes all scalar units; both GL
and SL retain -I. These are presentations of GL and SL, not PGL or PSL.

 D    Local ridges  Local ridge orbits  Oriented loops   Length:frequency
 5       75,15          11,4                180         3:82, 4:72, 6:26
 8       15,76           4,11               182         3:92, 4:14, 5:76

Each ridge is processed in both directions, without using local reduction
to skip coverage. All ledger entries are complete, with exact reverse-frame
witnesses. Global two-cell orbit minimization is not implemented; the complete
redundant set supplies the required attaching relations. Rank-four ridges all
have K-rank-two supporting vectors and explicit positive interior witnesses.

 D       GL raw       GL simplified       SL raw       SL simplified
 5      20 / 390          3 / 10          77 / 1564          2 / 6
 8      17 / 291          3 / 11          65 / 1168          2 / 6

Entries are numbers of generators / relators. Unsimplified presentations
retain finite stabilizer multiplication tables, all edge inclusions,
coherent tree lifts, exact matrix images and relator provenance. Both old/new
generator maps are retained for each Tietze simplification.

Edges are subdivided to handle inversions. The full setwise edge group
becomes the midpoint vertex group. D5 has no inversions. D8 has one inverted
edge, with oriented stabilizer order two and setwise stabilizer order four:
its endpoint-reversing matrices square to I. The algorithm checks their
actual squares in the oriented stabilizer and never infers reversal orders
from the order of the setwise group. A separate order-four regression checks
that a matrix whose square is -I retains that square in its multiplication
table.

Constructive matrix-to-word reduction crosses separating facets while
strictly decreasing an exact positive trace objective in a discrete rational
value set. It provides the scalar word for epsilon*I, where epsilon=w (D5)
or 1+w (D8), verified fundamental by Hecke without GRH. Adding this central
relation to GL precedes the index-four determinant-mod-squares Schreier
rewrite. Every generator retains its GL lift and its determinant-one scalar
correction. This avoids treating SL as finite index in the original GL.

Read results/TWO_CELL_RESULTS.txt for counts and timings, and for each D:

  results/two_cells/D5/FIRST_TWO_CELL.txt   first inspectable ordered cycle
  results/two_cells/D5/first_two_cell.json exact first-cycle certificate
  results/two_cells/D5/coverage.json       all ridges and both orientations
  results/two_cells/D5/loops/              exact replayable loop certificates
  results/two_cells/D5/GL_presentation.json
  results/two_cells/D5/SL_presentation.json
  results/two_cells/D5/GL_simplified.json   matrices and both Tietze maps
  results/two_cells/D5/SL_simplified.json
  results/two_cells/D5/PRESENTATIONS.txt   readable generators and relators

The same paths exist under D8. results/two_cell_summary.json records stage
timings and completion statuses. The saved baseline initial-coverage ledgers
in results/two_cells/baseline/ retain uncached loop timings; repeated runs
report the cost of exact cache replay instead. Timings within Julia include
first-use compilation and exclude process startup/package import. They are
not filesystem-cold benchmarks. No new backend/sanitizer audit was required.

PRESENTATION_JUSTIFICATION.md connects the complete arithmetic enumeration,
positive ridge cycles, full stabilizers, inversion handling and coherent
attaching paths to the contractible well-rounded complex and the group-action
presentation theorem. This supplies faithfulness as well as the exactly
checked homomorphism to matrices. Raw matrix-relator checks alone would not.

Licensing/attribution: see backend/NOTICE.txt, backend/UPSTREAM_COPYING and
backend/PERMUTALIB_COPYING.LGPL. No upstream sources have been relicensed.

Following the acceptance audit, the user requested three upstream reports.
They were submitted as polyhedral_common issues #13, #14 and #15. Exact bodies,
tested reproductions and verified publication receipts are in
upstream_reports/README.txt and upstream_reports/submission_manifest.json.

Field extension, constructive words and dataset
----------------------------------------------

The fields/words/data handoff extends the same engine to D12, D13 and D40.
New stage records, presentations and certificates are in results/fields/D*/.
The original D5/D8 outputs, executable and configuration remain in place;
results/baselines/two_cells_2026_10_17/ additionally preserves their exact
source/data archive and 454-file manifest. Backend changes for this extension
only admit the additional three input models; no enumeration algorithm was
changed. The new executable uses backend/config_fields.json and is selected
for D12/D13/D40. Reproduce its build/configuration with:

  python3 experiments/arithmetic_voronoi/scripts/build_backend.py \
    --target AV_ExportPerfectOrbits --configure-fields

  julia --startup-file=no --project=experiments/arithmetic_voronoi \
    experiments/arithmetic_voronoi/scripts/run_fields.jl 12 13 40

run_fields.jl finishes GL, SL and the separate PSL quotient. Stage statuses
and elapsed times are saved independently. A complete coverage cache is
replayed without repeatedly rewriting the large ledger. Failed or pending
stages cannot promote the raw presentation status to complete. Initial
complete run timings are retained separately from later cache replays.

CONSTRUCTIVE_WORDS_STATUS.md documents the public API and its guarantees;
KW_COMPARISON.md documents the four-field literature maps and the D8
transcription discrepancy. DATASET_SCHEMA.md describes exact encodings,
stable group identities, maps, certificates and the two verification levels.
DATASET_PROPOSAL.md is a local proposal, with no LMFDB submission or contact.

The full reproduction/verification command (with installed dependencies) is:

  julia --startup-file=no --project=experiments/arithmetic_voronoi \
    experiments/arithmetic_voronoi/scripts/reproduce_fields_words_data.jl

The extended test/runtests.jl includes all baseline tests plus the five-field
word, central-sign, membership, resource-exhaustion, cache-corruption,
Kirchheimer-Wolfart and nonprincipal-ideal regressions. Targeted sanitizer
results for the new input guard are kept separately in results/sanitizer_fields.
Their reproduction commands are:

  python3 experiments/arithmetic_voronoi/scripts/build_backend.py \
    --target AV_ExportPerfectOrbits --sanitizer
  python3 experiments/arithmetic_voronoi/scripts/run_sanitizer.py --fields
  julia --startup-file=no --project=experiments/arithmetic_voronoi \
    experiments/arithmetic_voronoi/scripts/validate_sanitizer.jl --fields

The dataset contains a copy of the pinned project, source and checker and can
be validated without an exporter executable. See DATASET_SCHEMA.md for the
fresh-process commands, including the stronger --complete check.
