ArithmeticVoronoi
=================

Exact arithmetic Voronoi computations through OSCAR: perfect forms,
GL(2,O), SL(2,O) and PSL(2,O) presentations, and verified matrix-to-word
reduction. Supported real quadratic field discriminants are 5, 8, 12, 13
and 40. In each case O is the maximal order and the lattice is O^2.
GL and SL retain -I; PSL is the separate quotient SL/{I,-I}.

After cloning, run the commands from the repository root, beside Project.toml.
No previous workspace, backend binary, configuration or computation cache is
required. Magma is not a dependency.

Requirements
-------------

For the saved data: Git, Git LFS, Python 3.12 or newer, and tar.
For verification and computations: Julia 1.13 (tested with 1.13.2).
Project.toml and Manifest.toml pin the Julia packages, including OSCAR 1.8.2
and Hecke 0.39.22. The instantiate commands below download missing packages
into your Julia depot without changing your default Julia environment.

For new perfect-form enumeration, the external backend additionally needs
Linux, a C++20 compiler, CMake, Ninja, Make, pkg-config, prlimit (from util-linux),
and development files for GMP, Eigen and Boost.Serialization. Install these
with your distribution's package manager. Initial setup needs network access.

Get the checkout and data
-------------------------

  git clone https://github.com/dimpase/arithmetic_voronoi.git
  cd arithmetic_voronoi
  git lfs install --local
  git lfs pull

If you already have a checkout, start with the last two commands in its root.
The two LFS archives hold the five-field dataset and D5/D8 regression fixtures.
They contain no upstream source snapshot or compiled backend.

Read and verify the saved results
---------------------------------

No C++ build is needed for this section. Verify and unpack the dataset:

  python3 -c 'import hashlib,json,pathlib; p=pathlib.Path("results/dataset"); m=json.loads((p/"archive_manifest.json").read_text()); assert hashlib.sha256((p/m["archive"]).read_bytes()).hexdigest()==m["sha256"]; print("archive checksum OK")'
  tar -xzf results/dataset/arithmetic_groups_dataset_v1.tar.gz -C results/dataset

Open results/dataset/arithmetic_groups_dataset_v1/fields/D12/PRESENTATIONS.txt
for a readable presentation. Replace D12 by D5, D8, D13 or D40 for another
field. The adjacent JSON files contain exact matrices, relations and maps.
results/FIELDS_WORDS_DATA_REPORT.txt summarizes all five completed cases.

Install the bundled checker's dependencies and verify the data:

  julia --startup-file=no --project=results/dataset/arithmetic_groups_dataset_v1/software \
    -e 'using Pkg; Pkg.instantiate(;allow_autoprecomp=false)'
  julia --startup-file=no --project=results/dataset/arithmetic_groups_dataset_v1/software \
    results/dataset/arithmetic_groups_dataset_v1/software/scripts/verify_dataset.jl \
    results/dataset/arithmetic_groups_dataset_v1

Add --complete to the last command to replay the geometric coverage and fresh
Tietze isomorphisms too. The default check verifies exact arithmetic, matrix
relations, maps, invariants and word certificates. Both modes work without
a backend executable. Keep the extracted bundle unchanged for its hash checks.

Install for new computations
-----------------------------

Install this checkout's pinned Julia dependencies:

  julia --startup-file=no --project=. \
    -e 'using Pkg; Pkg.instantiate(;allow_autoprecomp=false)'

Choose a writable cache outside the checkout, then configure and build the
external backend. Run this block in one shell; you may set AV_BACKEND_DIR to
another absolute path. Use the same value for all three Python commands.

  AV_BACKEND_DIR="${XDG_CACHE_HOME:-$HOME/.cache}/arithmetic-voronoi/backend"
  python3 scripts/configure_backend.py \
    --source "$AV_BACKEND_DIR/source" --build "$AV_BACKEND_DIR/build"
  python3 scripts/build_backend.py --target AV_ExportPerfectOrbits \
    --source "$AV_BACKEND_DIR/source" --build "$AV_BACKEND_DIR/build" --configure-default
  python3 scripts/build_backend.py --target AV_ExportPerfectOrbits \
    --source "$AV_BACKEND_DIR/source" --build "$AV_BACKEND_DIR/build" --configure-fields

The configure script fetches the exact upstream and submodule revisions from
backend/pins.json and applies our CMake integration patch. The build uses one
job and a 10 GiB virtual-memory ceiling. The second install command reuses
the compiled target and configures the additional fields.

These commands create backend/bin/ and backend/config*.json for YOUR checkout.
Both configurations are needed to run all five fields and the full test suite.
They contain local executable paths and hashes and are ignored by Git.
Reconfigure and rebuild after relocating the checkout; do not copy another
installation's configuration or executable. Logs are written to results/build/.
Upstream sources remain in the external cache, outside this repository.

Run a computation
------------------

Compute perfect forms, complete ridge coverage and GL/SL/PSL presentations:

  julia --startup-file=no --project=. scripts/run_fields.jl 12

Use multiple arguments, e.g. 12 13 40; no arguments selects those three fields.
D is a field discriminant: D8 means Q(sqrt(2)), D12 means Q(sqrt(3)), and
D40 means Q(sqrt(10)). Stages, timings, presentations and certificates appear
under results/fields/D12/. D12 has three perfect-form orbits and 228 oriented
ridge loops; its simplified SL presentation has two generators and six relators.

For a complete matrix-to-word example, including exact certificate replay:

  julia --startup-file=no --project=. scripts/word_example.jl 12

This uses its own results/word_example/D12/ directory and converts
[1 w; 0 1] to a word in the simplified SL generators. Words are signed,
one-based generator indices, multiplied in written order.

For D5/D8 runs that reuse the original regression geometry, restore it first:

  python3 scripts/restore_fixtures.py
  julia --startup-file=no --project=. scripts/run_fields.jl 5 8

The restoration script verifies all file hashes and refuses to overwrite
different existing files. The word example instead computes its own ridge
data and can be run for any of the five fields without restoring these fixtures.

Tests and reproducibility
--------------------------

After installing both backend configurations, run:

  python3 scripts/restore_fixtures.py
  julia --startup-file=no --project=. test/runtests.jl

The original delivery passed 3806 assertions; its receipts and dataset remain
versioned. Fresh runs create ignored caches and scratch data and may update
small tracked summary reports. First use includes Julia compilation. See
run.json for per-stage timings and cache use; saved warm timings are not
cold end-to-end benchmarks.

If an archive is a short Git LFS pointer file, run git lfs pull. A
missing_executable error means the relevant backend configuration has not
been built here. Changed executables invalidate their caches. For a separate
D8 computation, use run_field_pipeline(8; directory="results/my_D8",
reuse_baseline=false) after using ArithmeticVoronoi in Julia. See
USER_GUIDE.txt for the full API and resource-limit behaviour.

Further reading
----------------

  USER_GUIDE.txt                 API examples, output files and troubleshooting
  PRESENTATION_JUSTIFICATION.md  completeness and group-action conventions
  CONSTRUCTIVE_WORDS_STATUS.md   word reduction, termination and certificates
  DATASET_SCHEMA.md             exact encodings and verification levels
  KW_COMPARISON.md               explicit Kirchheimer-Wolfart comparisons
  BACKEND_AUDIT.md               historical engineering acceptance record
  backend/NOTICE.txt            external dependency and attribution details
  LICENSE                       project license

The historical audit and saved receipts describe their original execution
environment; use this README and USER_GUIDE.txt for installation commands.
