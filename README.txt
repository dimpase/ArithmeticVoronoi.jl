ArithmeticVoronoi
=================

Compute presentations and matrix-to-word reductions for GL(2,O), SL(2,O)
and PSL(2,O) in Julia, using OSCAR and either the native Julia perfect-form
enumerator or the PerfectForms C++ backend. O is the maximal order of a
quadratic field, acting on the free lattice O^2. The native backend also
handles imaginary quadratic fields, using Hermitian forms.
docs/native_survey.md records the 16 real fields; docs/imaginary_quadratic.md
gives five imaginary examples, results and comparisons with the literature.
GL and SL retain -I; PSL is the separate quotient SL/{I,-I}.

Start a new group computation
-----------------------------

These instructions compute presentations from perfect forms and ridge cycles.
Saved datasets, Git LFS and Python are optional; they are used only in the
archive and regression-test sections below. Julia 1.13 is required (tested
with 1.13.2). Magma is not a dependency.

1. Enter the checkout and install its Julia dependencies

For a new checkout, skip downloading the optional saved-data archives:

  GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/dimpase/arithmetic_voronoi.git
  cd arithmetic_voronoi

For an existing checkout, cd to its root instead. Run all following shell
commands beside Project.toml:

  julia --startup-file=no --project=. \
    -e 'using Pkg; Pkg.instantiate(;allow_autoprecomp=false)'

Project.toml and Manifest.toml select the tested OSCAR, Hecke and PerfectForms
versions. This installs into your Julia depot and leaves your default Julia
environment unchanged. The first import and first computation compile Julia
code and can pause for several minutes before any stage output appears.

2. Choose the perfect-form enumerator

The examples below use backend=:native. This computes a starting perfect form
and the complete neighbour graph using Hecke and OSCAR, directly in Julia.
No C++ exporter build or PERFECT_FORMS_BACKEND_CONFIG setting is needed.
See docs/native_backend.md for the algorithm, supported scope and comparison
tests. The existing C++ path remains the API default for compatibility with
saved computations; pass backend=:native explicitly to select the new path.

To use C++ instead, select an existing PerfectForms installation, for example:

  export PERFECT_FORMS_BACKEND_CONFIG="$HOME/work/tmp/PerfectForms.jl/backend/config.json"

It must be the absolute path to the installed config.json FILE. Follow "Build
the optional C++ backend" below if needed, and use backend=:polyhedral_common
in the Julia call. No separate Oscar, Hecke or GAP installation is needed beyond
Pkg.instantiate() for either enumerator.

3. Start Julia and compute the presentations

In the shell:

  julia --startup-file=no --project=.

At the Julia prompt:

  using ArithmeticVoronoi, Oscar

  R = run_field_pipeline(12; backend=:native, progress=true,
      directory="results/my_native_D12", reuse_baseline=false);

This computes GL(2,O), SL(2,O) and PSL(2,O) for O = Z[sqrt(3)]. It enumerates
and validates perfect forms, computes the two-cells, constructs all three
presentations and simplifies them with certified generator maps. Stage names
are printed as it runs. The trailing semicolon avoids printing the large
result object.

For the optional C++ path, if Julia was already running when you set the shell
variable, set it in that Julia session too, before the computation:

  ENV["PERFECT_FORMS_BACKEND_CONFIG"] =
      expanduser("~/work/tmp/PerfectForms.jl/backend/config.json")

D is the field discriminant: 5, 8, 12, 13 and 40 correspond to Q(sqrt(5)),
Q(sqrt(2)), Q(sqrt(3)), Q(sqrt(13)) and Q(sqrt(10)). Change 12 and the output
directory together for a different field. These are discriminants, not always
radicands: for example use 60 for Q(sqrt(15)).

For a Bianchi group, use a negative fundamental discriminant and the native
backend, for example run_field_pipeline(-7; backend=:native,
directory="results/my_Dminus7", reuse_baseline=false). Use -4 for Q(i), not -1.
See docs/imaginary_quadratic.md for the complete examples and Page comparison.

reuse_baseline=false computes ridge data in your run directory without using
the archived D5/D8 fixtures. Completed compatible data in that directory and
validated C++ perfect-form caches can still be reused on subsequent runs.
Native perfect forms are computed afresh. Use a separate output directory for
each backend: their representative ordering and generator labels can differ.

To compute only perfect forms and their labelled neighbours:

  P = arithmetic_problem(13);
  E = enumerate_perfect_forms(P; backend=:native, progress=true);
  length(E.representatives)  # 9
  [r.stabilizer_order for r in E.representatives]
  E.neighbours

Work with a group in OSCAR
--------------------------

Continue in the same Julia session. Select the simplified SL presentation,
inspect its exact matrices and relators, and construct the OSCAR group:

  S = R.simplified[:SL];
  S.generator_matrices
  S.relators

  F = free_group(length(S.generator_matrices));
  G, toG = quo(F, [F(r) for r in S.relators]);
  gens(G)
  relators(G)

G is now an OSCAR FPGroup representing SL(2,O). Generator i corresponds to
S.generator_matrices[i]. Each relator is a vector of signed, one-based generator
indices: [1,-2,1] means g1 * inv(g2) * g1 and [] is the identity.

For example, compute the abelianization and its quotient map:

  A, toA = maximal_abelian_quotient(G);
  elementary_divisors(abelian_invariants(A))   # D12: [2, 6], so C2 x C6
  order(A)                # 12

G itself is infinite. A is its finite abelian quotient in this example.
To work with GL or PSL instead, select R.simplified[:GL] or R.simplified[:PSL]
and repeat the construction of F and G. The PSL matrices are determinant-one
lifts: its relators hold modulo +/-I. The raw presentations are available as
R.gl, R.sl and R.psl, with generator_matrices and relators fields.

Convert a matrix into a group element
-------------------------------------

Using the SL group G, toG and S constructed above:

  P = R.enumeration.problem;
  g = matrix(P.field, [1 P.generator; 0 1]);
  c = matrix_to_word(R, g; target=:SL);
  word = Int.(c["word"])
  x = toG(F(word))

x is the element of G represented by the input matrix g. P.generator is the
integral field generator w (sqrt(3) for D12). Evaluate its image in the
abelianization, or replay the exact matrix and reduction certificate:

  toA(x)
  map_word(x, S.generator_matrices; init=identity_matrix(P.field, 2)) == g
  verify_word_certificate(R, g, c)   # true

The word and generator matrices must come from the same run and alphabet.
The default alphabet is simplified; set simplified=false for a raw one.
target=:GL accepts integral matrices with unit determinant. target=:SL and
:PSL require determinant-one matrices; PSL represents the matrix modulo +/-I.
The returned words are verified, but need not be shortest or canonical.

Run from the shell
-------------------

Run the native group pipeline directly from the shell:

  julia --startup-file=no --project=. -e 'using ArithmeticVoronoi; run_field_pipeline(12; backend=:native, progress=true, directory="results/my_native_D12", reuse_baseline=false);'

For the full 16-field native batch (5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 37,
40, 41, 44, 56, 60), including GL/SL/PSL and matrix-to-word checks:

  julia --threads=1 --startup-file=no --project=. scripts/run_native_fields.jl

This writes results/native_fields/D*/ and an incrementally updated summary.json.
Add --resume to retain fields already marked complete in that summary.
Select fields and a different output directory, for example:

  julia --threads=1 --startup-file=no --project=. scripts/run_native_fields.jl \
    --output=results/my_native_batch 17 21 60

The older convenience scripts use the default C++ backend and its
environment variable:

  julia --startup-file=no --project=. scripts/run_fields.jl 12
  julia --startup-file=no --project=. scripts/word_example.jl 12

run_fields.jl writes to results/fields/D12/; multiple arguments such as
5 8 12 13 are accepted, and no arguments selects 12 13 40. Archived fixtures
are optional even for D5/D8. They are reused only if available and compatible.
word_example.jl uses results/word_example/D12/, computes its own ridge data,
converts [1 w; 0 1] to an SL word, and replays its certificate.

Each pipeline run writes field.json, enumeration.json, GL_raw.json,
SL_raw.json, PSL_raw.json, the corresponding *_simplified.json files, and
run.json with stages, timings and cache use. A separate run with
reuse_baseline=false stores its geometric certificates under two_cells/.
In Julia, inspect R.run["presentations"] for the presentation sizes and
abelianizations, and R.run["geometry"] for the geometric counts.

Build the optional C++ backend
------------------------------

PerfectForms provides the exporter and its CMake build. This needs Linux or
macOS, Git, a C++20 compiler and standard library with std::format, CMake 3.24+,
Make or Ninja, pkg-config, and development files for GMP, Eigen 3,
Boost.Serialization and nauty/Traces. Linux runtime limits use prlimit from
util-linux; macOS uses the system ps command. Initial setup needs network
access. See https://github.com/dimpase/PerfectForms.jl for dependency details.

On macOS install current Xcode Command Line Tools and run:

  brew install cmake ninja pkgconf gmp eigen@3 boost nauty
  export CMAKE_PREFIX_PATH="$(brew --prefix eigen@3):$(brew --prefix)"

After the Julia dependency installation in step 1, run from this repository's
root. These commands also work when Julia's package source is read-only:

  perfect_forms_source="$(julia --startup-file=no --project=. -e 'using PerfectForms; print(pkgdir(PerfectForms))')"
  cmake -S "$perfect_forms_source" -B build/perfect-forms \
    -DCMAKE_BUILD_TYPE=Release '-DCMAKE_CXX_FLAGS_RELEASE=-O1 -DNDEBUG' \
    -DCMAKE_INSTALL_PREFIX="$PWD/backend-install"
  cmake --build build/perfect-forms --parallel 1
  cmake --install build/perfect-forms --prefix "$PWD/backend-install"
  export PERFECT_FORMS_BACKEND_CONFIG="$PWD/backend-install/config.json"

CMake fetches and verifies the pinned upstream sources and recursive submodules.
The first configure and build can take several minutes. Build serially to limit
compiler memory; on Linux, optionally prefix the build command with
prlimit --as=10737418240 for a hard 10 GiB address-space ceiling.

One exporter supports all these fields. config.json records its absolute
executable path and SHA256; backend-build.json records build provenance.
Reinstall after rebuilding or moving the installation. Set the configuration
variable in each new shell, or in ENV inside the Julia session.

For perfect-form enumeration in other dimensions or arbitrary T-spaces, use
PerfectForms directly. ArithmeticVoronoi's native group pipeline handles rank
two over maximal real quadratic orders, computing units and class groups with
Hecke rather than a table of field-specific choices.

Optional saved results and regression fixtures
-----------------------------------------------

results/native_presentations.json contains all 48 simplified GL/SL/PSL matrix
presentations from the 16-field native survey. See docs/native_survey.md for
counts, timings and Julia commands to load them directly into OSCAR.

Download the two frozen data archives with Git LFS when you want to inspect
or verify that delivery:

  git lfs install --local
  git lfs pull

results/FIELDS_WORDS_DATA_REPORT.txt summarizes the saved five-field dataset.
USER_GUIDE.txt's archive section gives checksum, extraction and verification
commands; no backend is required to verify those saved results.

The native Julia tests need no exporter or archived fixtures:

  julia --startup-file=no --project=. test/native.jl
  julia --startup-file=no --project=. test/quadratic_units.jl

To compare freshly computed native forms and labelled neighbours with C++,
set PERFECT_FORMS_BACKEND_CONFIG and run:

  julia --startup-file=no --project=. test/native_comparison.jl

The full regression suite also uses the original D5/D8 fixtures. With Python
3.12+ and the LFS archives available, restore the fixtures, then run the tests:

  python3 scripts/restore_fixtures.py
  julia --startup-file=no --project=. test/runtests.jl

The restoration script checks hashes and refuses to overwrite different
existing files. The C++ exporter's own tests can also be run with:

  ctest --test-dir build/perfect-forms --output-on-failure

First-use timings include Julia compilation. Saved warm timings are not cold
end-to-end benchmarks. Tests create ignored scratch data and may update small
tracked summary reports. CI builds and tests on Linux and macOS:
https://github.com/dimpase/arithmetic_voronoi/actions/workflows/ci.yml

Further reading
----------------

  docs/native_backend.md        native enumeration, scope and comparison tests
  docs/native_survey.md         16-field counts, timings and full presentations
  USER_GUIDE.txt                 API details, archives and troubleshooting
  PRESENTATION_JUSTIFICATION.md  completeness and group-action conventions
  CONSTRUCTIVE_WORDS_STATUS.md   word reduction, termination and certificates
  DATASET_SCHEMA.md             exact encodings and verification levels
  KW_COMPARISON.md               explicit Kirchheimer-Wolfart comparisons
  BACKEND_AUDIT.md               historical engineering acceptance record
  DEPENDENCIES.txt               backend package and attribution details
  LICENSE                       project license

Historical receipts describe their original execution environment; use this
README and USER_GUIDE.txt for current installation and computation commands.
