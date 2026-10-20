ArithmeticVoronoi
=================

Compute presentations and matrix-to-word reductions for GL(2,O), SL(2,O)
and PSL(2,O) in Julia, using OSCAR and the PerfectForms C++ backend. O is the
maximal order of a real quadratic field, acting on the free lattice O^2.
The complete group pipeline supports field discriminants 5, 8, 12, 13 and 40.
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

2. Select the installed PerfectForms backend

If you have already built and installed PerfectForms at
~/work/tmp/PerfectForms.jl/backend, reuse that installation:

  export PERFECT_FORMS_BACKEND_CONFIG="$HOME/work/tmp/PerfectForms.jl/backend/config.json"

Adjust the path for your installation. It must be the absolute path to the
installed config.json FILE. If you need a new backend, follow "Build a backend"
below, then return here. No separate Oscar, Hecke or GAP installation is needed
beyond Pkg.instantiate().

3. Start Julia and compute the presentations

In the shell:

  julia --startup-file=no --project=.

At the Julia prompt:

  using ArithmeticVoronoi, Oscar
  default_backend().executable

  R = run_field_pipeline(12; directory="results/my_D12", reuse_baseline=false);

This computes GL(2,O), SL(2,O) and PSL(2,O) for O = Z[sqrt(3)]. It enumerates
and validates perfect forms, computes the two-cells, constructs all three
presentations and simplifies them with certified generator maps. Stage names
are printed as it runs. The trailing semicolon avoids printing the large
result object.

If Julia was already running when you set the shell variable, set it in that
Julia session too, before the computation:

  ENV["PERFECT_FORMS_BACKEND_CONFIG"] =
      expanduser("~/work/tmp/PerfectForms.jl/backend/config.json")

D is the field discriminant: 5, 8, 12, 13 and 40 correspond to Q(sqrt(5)),
Q(sqrt(2)), Q(sqrt(3)), Q(sqrt(13)) and Q(sqrt(10)). Change 12 and the output
directory together for a different field. D40 is a substantially larger run.

reuse_baseline=false computes ridge data in your run directory without using
the archived D5/D8 fixtures. Completed compatible data in that directory and
validated perfect-form caches can still be reused on subsequent runs.

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

  A, toA = maximal_abelian_quotient(FinGenAbGroup, G);
  elementary_divisors(A)   # D12: [2, 6], so C2 x C6
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

With the same backend environment variable set, these scripts run complete
computations without an interactive session:

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

Build a backend (if you do not already have one)
------------------------------------------------

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

One exporter supports all five fields. config.json records its absolute
executable path and SHA256; backend-build.json records build provenance.
Reinstall after rebuilding or moving the installation. Set the configuration
variable in each new shell, or in ENV inside the Julia session.

For perfect-form enumeration in other fields, dimensions or T-spaces, use
PerfectForms directly. ArithmeticVoronoi's complete group pipeline currently
supports the five discriminants listed above.

Optional saved results and regression fixtures
-----------------------------------------------

Download the two frozen data archives with Git LFS when you want to inspect
or verify that delivery:

  git lfs install --local
  git lfs pull

results/FIELDS_WORDS_DATA_REPORT.txt summarizes the saved five-field dataset.
USER_GUIDE.txt's archive section gives checksum, extraction and verification
commands; no backend is required to verify those saved results.

The regression suite uses the original D5/D8 fixtures. With Python 3.12+ and
the LFS archives available, restore the fixtures, then run the tests:

  python3 scripts/restore_fixtures.py
  julia --startup-file=no --project=. test/runtests.jl

The restoration script checks hashes and refuses to overwrite different
existing files. Native backend tests can also be run with:

  ctest --test-dir build/perfect-forms --output-on-failure

First-use timings include Julia compilation. Saved warm timings are not cold
end-to-end benchmarks. Tests create ignored scratch data and may update small
tracked summary reports. CI builds and tests on Linux and macOS:
https://github.com/dimpase/arithmetic_voronoi/actions/workflows/ci.yml

Further reading
----------------

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
