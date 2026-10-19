ArithmeticVoronoi
=================

Exact arithmetic Voronoi computations through OSCAR: perfect forms,
GL(2,O), SL(2,O) and PSL(2,O) presentations, and verified matrix-to-word
reduction. The complete presentation and word pipeline currently supports
real quadratic field discriminants 5, 8, 12, 13 and 40. In each case O is
the maximal order and the lattice is O^2.
GL and SL retain -I; PSL is the separate quotient SL/{I,-I}.

For perfect-form enumeration in other fields, dimensions or T-spaces, use
PerfectForms directly: https://github.com/dimpase/PerfectForms.jl
Its interface passes general upstream T-space descriptions to the enumerator.

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

The exporter is provided by the separate PerfectForms Julia package:
https://github.com/dimpase/PerfectForms.jl
Its README gives dependency installation and CMake options, including macOS
Homebrew prefixes and generic system-nauty detection.

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

Locate the pinned PerfectForms dependency and build it with ordinary CMake:

  perfect_forms_source="$(julia --startup-file=no --project=. -e 'using PerfectForms; print(pkgdir(PerfectForms))')"
  cmake -S "$perfect_forms_source" -B build/perfect-forms -G Ninja \
    -DCMAKE_BUILD_TYPE=Release '-DCMAKE_CXX_FLAGS_RELEASE=-O1 -DNDEBUG' \
    -DCMAKE_INSTALL_PREFIX="$PWD/backend-install"
  cmake --build build/perfect-forms --parallel 1
  ctest --test-dir build/perfect-forms --output-on-failure
  cmake --install build/perfect-forms
  export PERFECT_FORMS_BACKEND_CONFIG="$PWD/backend-install/config.json"

The build directory and install prefix can be any writable directories; the
Julia dependency's source directory need not be writable. CMake fetches and
verifies PerfectForms' exact upstream and submodule pins. An installed
nauty/Traces is required, using pkg-config metadata and standard CMake
searches. No source tree or binary is bundled in ArithmeticVoronoi.

On macOS install current Xcode Command Line Tools, then, with Homebrew:

  brew install cmake ninja pkgconf gmp eigen@3 boost nauty git-lfs python

Add -DCMAKE_PREFIX_PATH="$(brew --prefix eigen@3);$(brew --prefix)" to the
configure command above. These prefixes work on Apple Silicon and Intel Macs.

Build serially to limit compiler memory. On Linux, prefix the build command
with prlimit --as=10737418240 for a hard 10 GiB address-space ceiling.
One exporter supports all five fields. backend-install/config.json selects
its absolute path and SHA256; backend-install/backend-build.json records
its compiler, upstream pins and system-nauty provenance. Set the environment
variable in each new shell, or pass BackendConfig with an absolute executable
path. Reinstall after changing or moving a binary; its SHA controls cache use.

PerfectForms owns our C++ code, CMake build, raw JSON protocol and process
controls. ArithmeticVoronoi owns the arithmetic models, exact validation,
mathematical caches, two-cells, presentations and constructive words.

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

  ctest --test-dir build/perfect-forms --output-on-failure
  python3 scripts/restore_fixtures.py
  julia --startup-file=no --project=. test/runtests.jl

The original delivery passed 3806 assertions; its receipts and dataset remain
versioned. Fresh runs create ignored caches and scratch data and may update
small tracked summary reports. First use includes Julia compilation. See
run.json for per-stage timings and cache use; saved warm timings are not
cold end-to-end benchmarks.

.github/workflows/ci.yml builds the exporter and runs native and Julia tests
on Linux and macOS. Current logs and platform results are available at:
https://github.com/dimpase/arithmetic_voronoi/actions/workflows/ci.yml

If an archive is a short Git LFS pointer file, run git lfs pull. A
missing_executable error means the backend has not been installed here;
run cmake --install build/perfect-forms after a successful build, and set
PERFECT_FORMS_BACKEND_CONFIG as above. For a separate
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
  DEPENDENCIES.txt               backend package and attribution details
  LICENSE                       project license

The historical audit and saved receipts describe their original execution
environment; use this README and USER_GUIDE.txt for installation commands.
