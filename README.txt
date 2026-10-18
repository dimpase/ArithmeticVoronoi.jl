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
Linux or macOS, a C++20 compiler and standard library (including std::format),
CMake 3.24 or newer, Ninja, pkg-config, and development
files for GMP, Eigen, Boost.Serialization and nauty/Traces. Install nauty's
development headers and library before configuring the backend.
Install these with your distribution's package manager. Initial setup needs
network access. Python is not required to build or install the exporter.
Linux computations also use prlimit from util-linux. macOS computations use
the system ps command to monitor resident memory.

On macOS, install current Xcode Command Line Tools and the Homebrew packages:

  xcode-select --install
  brew install cmake ninja pkgconf gmp eigen@3 boost nauty git-lfs python

Configure with Homebrew's prefixes so CMake finds Eigen 3, including its
versioned installation, on both Apple Silicon and Intel Macs:

  cmake --preset release \
    -DCMAKE_PREFIX_PATH="$(brew --prefix eigen@3);$(brew --prefix)"

Continue with the build/install commands below. CMake reports a clear error
if the selected C++ standard library lacks std::format. See the package lists:
https://formulae.brew.sh/formula/nauty and https://formulae.brew.sh/formula/eigen@3

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

Configure, build and install the native exporter with CMake:

  cmake --preset release
  cmake --build --preset release
  cmake --install build

The release preset uses Ninja, one build job and -O1 to limit compiler memory.
It installs into this checkout's backend/ directory. For a hard 10 GiB build
memory ceiling on Linux, use prlimit --as=10737418240 before the build command.
No administrator access is needed. The equivalent configuration without presets is:

  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    '-DCMAKE_CXX_FLAGS_RELEASE=-O1 -DNDEBUG' -DCMAKE_INSTALL_PREFIX="$PWD/backend"
  cmake --build build --parallel 1
  cmake --install build

CMake fetches and verifies the exact upstream and submodule revisions in
backend/pins.json. Sources live in build/_deps/ and remain unmodified; the
upstream standalone build is not configured. No dependency source is tracked
in this repository. The native build and all generated files are ignored.

CMake requires an installed nauty/Traces and checks that its Traces interface
compiles and links. Configuration stops with an installation hint if it is
missing or unusable. For a nonstandard installation, set CMAKE_PREFIX_PATH,
or set Nauty_INCLUDE_DIR (traces.h), Nauty_CORE_INCLUDE_DIR (nauty.h), and
Nauty_LIBRARY explicitly. Detection uses pkg-config metadata when present
and standard CMake header/library searches. Stock Traces uses the
backend's exact group construction from generators.

One install creates backend/bin/AV_ExportPerfectOrbits, both backend/config.json
and backend/config_fields.json, and backend/backend-build.json. Both selectors
use the same executable, supporting all five fields. They record absolute
paths and executable hashes. Reinstall after rebuilding or relocating the
checkout; changed executables invalidate computation caches.

An external build directory works too: replace build in the explicit commands
with a writable directory of your choice. For an existing clean pinned source
checkout, set -DFETCHCONTENT_SOURCE_DIR_POLYHEDRAL_COMMON=/absolute/source/path
at configuration. Its commit and submodule revisions are still checked.

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
The pipeline reuses restored fixtures only when their exact source enumeration
matches; otherwise it computes and validates fresh two-cell certificates.

Tests and reproducibility
--------------------------

Check the native exporter, then run the Julia regressions:

  ctest --preset release
  python3 scripts/restore_fixtures.py
  julia --startup-file=no --project=. test/runtests.jl

The original delivery passed 3806 assertions; its receipts and dataset remain
versioned. Fresh runs create ignored caches and scratch data and may update
small tracked summary reports. First use includes Julia compilation. See
run.json for per-stage timings and cache use; saved warm timings are not
cold end-to-end benchmarks.

.github/workflows/ci.yml builds the exporter and runs native and Julia tests
on Linux and macOS. Local validation was performed on Linux; a macOS runner
is needed to establish the macOS build and runtime results.

If an archive is a short Git LFS pointer file, run git lfs pull. A
missing_executable error means the backend has not been installed here;
run cmake --install build after a successful build. For a separate
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
