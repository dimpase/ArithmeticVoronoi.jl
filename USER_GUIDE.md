# ArithmeticVoronoi user guide

ArithmeticVoronoi computes perfect forms, finite presentations and verified
matrix words for GL(2,O), SL(2,O) and PSL(2,O). Here O is a maximal quadratic
order and the lattice is the free module O². Both real and imaginary fields
are supported by the native Julia backend. Inputs are not restricted to the
fields already in our surveys.

The group algorithm follows Braun, Coulangeon, Nebe and Schönnenbeck,
[*Computing in arithmetic groups with Voronoi's algorithm*](https://arxiv.org/abs/1407.6234).
For enumeration details see the [native backend guide](docs/native_backend.md);
for the group-action argument see
[presentation justification](PRESENTATION_JUSTIFICATION.md).

## Install and start Julia

Use Julia 1.13; the tested installation uses Julia 1.13.2, OSCAR 1.8.2 and
Hecke 0.39.22. Run shell commands from the checkout root, beside `Project.toml`:

```sh
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/dimpase/arithmetic_voronoi.git
cd arithmetic_voronoi
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(;allow_autoprecomp=false)'
julia --startup-file=no --project=.
```

For an existing checkout, start at the `cd` step. The project and manifest
select the tested Julia dependencies; `--project=.` selects this environment.
The first import and first computation compile code and can take several
minutes before computation progress appears.

The native path needs no separately installed GAP, Magma, C++ compiler or
backend configuration. OSCAR brings GAP.jl and its GAP runtime. PerfectForms
remains a Julia dependency for the optional C++ interface; the pinned revision
does not build C++ during installation. Network access is needed to download
dependencies. Saved datasets and Git LFS are optional for new computations.

## Fields and group conventions

Pass a **fundamental field discriminant**, not a square-root radicand:

| Input D | Field | Integral generator w |
| ---: | --- | --- |
| 5 | Q(√5) | (1+√5)/2 |
| 8 | Q(√2) | √2 |
| 12 | Q(√3) | √3 |
| 40 | Q(√10) | √10 |
| 60 | Q(√15) | √15 |
| −3 | Q(√−3) | (1+√−3)/2 |
| −4 | Q(i) | i |
| −7 | Q(√−7) | (1+√−7)/2 |

The action is K-linear. GL includes all scalar units, and GL and SL retain
−I. PSL is the separate quotient SL/{I,−I}; its matrices are determinant-one
lifts. PGL, the quotient of GL by all scalar units, is a different group and
is not an output of this pipeline. Other lattices, nonmaximal orders and
higher ranks are outside the current native group interface.

## Compute presentations from scratch

At the Julia prompt:

```julia
using ArithmeticVoronoi, Oscar

R = run_field_pipeline(12; backend=:native, progress=true,
    directory="results/my_D12", reuse_baseline=false);
```

This enumerates perfect forms and neighbours, computes ridge cycles,
constructs GL/SL/PSL presentations, simplifies them with tracked generator
maps, and computes abelianizations. The trailing semicolon avoids printing
the large result object. `progress=true` reports enumeration and pipeline
stages. The total number of perfect-form orbits is discovered during traversal.

For an imaginary quadratic field:

```julia
Rim = run_field_pipeline(-7; backend=:native, progress=true,
    directory="results/my_Dminus7", reuse_baseline=false);
```

The API default remains `backend=:polyhedral_common` for compatibility with
older computations. Use `backend=:native` in new Julia-only runs, including
imaginary fields. The optional C++ comparison path currently covers the real
quadratic examples.

### Work with the group in OSCAR

Continue with the D12 result `R`:

```julia
S = R.simplified[:SL];
S.generator_matrices
S.relators

F = free_group(length(S.generator_matrices));
G, toG = quo(F, [F(r) for r in S.relators]);
gens(G)
relators(G)

A, toA = maximal_abelian_quotient(G);
elementary_divisors(abelian_invariants(A))   # [2, 6]: C2 × C6
order(A)                                   # 12
```

`G` is the infinite group SL(2,O); `A` is its abelianization, finite in this
example. Generator `i` corresponds to `S.generator_matrices[i]`. Use
`R.simplified[:GL]` or `R.simplified[:PSL]` and reconstruct `F` and `G` for
the other groups. The raw presentations are `R.gl`, `R.sl` and `R.psl`.

Use the one-argument `maximal_abelian_quotient(G)` shown above for both finite
and infinite abelianizations. In the returned invariant factors, `0` denotes
a copy of Z. For example, PSL over D = −7 has abelianization Z × C2.

### Express a matrix in these generators

Using the same `R`, `S`, `F`, `G`, `toG` and `toA`:

```julia
P = R.enumeration.problem;
g = matrix(P.field, [1 P.generator; 0 1]);
c = matrix_to_word(R, g; target=:SL);
word = Int.(c["word"]);
x = toG(F(word));
toA(x)

map_word(x, S.generator_matrices; init=identity_matrix(P.field, 2)) == g
verify_word_certificate(R, g, c)   # true
```

`P.generator` is the integral generator w from the field table. Words use
signed, one-based indices in written product order: `[1,-2,1]` means
g₁ g₂⁻¹ g₁, and `[]` is the identity. Words and matrices must belong to the
same run and generator alphabet.

`target=:GL` accepts integral matrices with unit determinant; `:SL` and
`:PSL` require determinant one. PSL certificates record equality modulo
±I and the exact lift sign. Words are verified, but need not be shortest or
canonical. Semilinear matrices and matrices outside the integral group are
rejected.

Use `simplified=false` for raw generators, `certificate=false` to return
only the word, or `cache_path="results/my_word.json"` to save the certificate.
The word reducer has optional `max_steps` and `timeout` budgets; these are
separate from native enumeration, which has no process deadline or memory
limit. See [word certificates](CONSTRUCTIVE_WORDS_STATUS.md) for details.

### Compute only perfect forms

```julia
P = arithmetic_problem(13);
E = enumerate_perfect_forms(P; backend=:native, progress=true);
length(E.representatives)                 # 9
E.completion_status                      # :complete
[r.stabilizer_order for r in E.representatives]
E.neighbours
```

Negative discriminants select Hermitian forms. `classical_problem()` gives
the binary rational case. The native enumerator uses exact lattice routines
in Hecke and polyhedral geometry in OSCAR; its matrices stay in Julia.

## Batches, output and restarting

The enumerator, presentation pipeline and HAP provider share a request record
describing the group, rank, field, order and lattice. `voronoi_request(D)` is
the shortcut for our quadratic GL experiments. See
[request records](docs/requests.md) for the full Julia/GAP schema and JSON examples.

Run a single field from the shell:

```sh
julia --startup-file=no --project=. -e 'using ArithmeticVoronoi; run_field_pipeline(12; backend=:native, progress=true, directory="results/my_D12", reuse_baseline=false);'
```

The real survey covers D = 5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 37, 40,
41, 44, 56, 60. The imaginary survey covers D = −3, −4, −7, −15, −23:

```sh
julia --threads=1 --startup-file=no --project=. scripts/run_native_fields.jl
julia --threads=1 --startup-file=no --project=. scripts/run_imaginary_fields.jl
```

Both batches include GL/SL/PSL presentations and exact matrix-to-word checks.
The imaginary script also loads HAP for comparison. Select a subset and
destination, or resume fields already marked complete:

```sh
julia --threads=1 --startup-file=no --project=. scripts/run_native_fields.jl --output=results/my_real_batch 17 21 60
julia --threads=1 --startup-file=no --project=. scripts/run_imaginary_fields.jl --output=results/my_imaginary_batch -4 -23
julia --threads=1 --startup-file=no --project=. scripts/run_native_fields.jl --resume
```

Both runners also accept `--requests=examples/voronoi_requests.json`; the
imaginary runner requires a file containing only imaginary GL requests.
The canonical request is saved with each result and checked when resuming.
To compute fresh complexes and higher homology through HAP's general interface:

```sh
julia --project=. scripts/run_hap.jl --requests=examples/voronoi_requests.json --degree=2 --output=results/my_homology
```

The default destinations are `results/native_fields/` and
`results/imaginary_fields/`. Each has an incrementally updated `summary.json`
and a `D*/` directory per field:

| File | Contents |
| --- | --- |
| `field.json` | Field, units, class group and action conventions |
| `enumeration.json` | Perfect forms, stabilizers, facets and neighbours |
| `two_cells/coverage.json` | Ridge coverage and references to loop certificates |
| `GL_raw.json`, `SL_raw.json`, `PSL_raw.json` | Raw presentations and geometric maps |
| `*_simplified.json` | Simplified matrices, relators and both generator maps |
| `run.json` | Stage status, timings, presentation sizes and cache use |
| `words/`, `word_validation.json` | Batch word certificates and their checks |

Inspect `R.run["presentations"]` and `R.run["geometry"]` interactively.
An ordinary `run_field_pipeline` call does not automatically run every word
test or literature comparison; the batch scripts add those steps.

`reuse_baseline=false` avoids the archived D5/D8 ridge fixtures. It does not
disable replay of compatible ridge and simplification certificates already
in the selected directory. Native perfect forms are enumerated afresh;
native enumeration checkpoint loading is not implemented. Use a **new output
directory** for a completely fresh computation and separate directories for
different backends: representative order and generator labels can differ.

The default `reuse_baseline=:auto` reuses archived fixtures only when their
exact source data match; `true` requires a matching baseline. Completed
certificates are validated on reuse. `--resume` acts at the batch level and
retains fields already marked complete in its summary.

## Load saved presentations

Ordinary Git includes compact bundles with 48 real-field presentations,
15 native imaginary-field presentations, and the 15 Page presentations:

- [Real survey bundle](results/native_presentations.json)
- [Imaginary survey bundle](results/imaginary_presentations.json)
- [Page comparison bundle](results/page_presentations.json)

For example, load the saved PSL group for D = −23 without recomputing geometry:

```julia
using Oscar, JSON3

data = JSON3.read(read("results/imaginary_presentations.json", String));
entry = data.fields["-23"].groups["PSL"];
F = free_group(length(entry.generator_matrices));
G, toG = quo(F, [F(Int.(r)) for r in entry.relators]);
A, toA = maximal_abelian_quotient(G);
elementary_divisors(abelian_invariants(A))   # [6, 0, 0, 0]
```

Matrix entries in the native bundles are pairs of rational strings `[a,b]`,
meaning a+bw with w=(s+√D)/2 and s=D mod 4. The Page bundle uses the same
matrix/word encoding but has a flat list of group records. The compact
bundles suffice for abstract group computations and exact matrix checks;
geometric matrix-to-word reduction needs a pipeline result `R` with its
geometry. See the survey guides for per-run file examples and timing scope.

## HAP, Magma and independent comparisons

HAP is a **GAP package**, available from Julia through **GAP.jl**. OSCAR
already includes GAP.jl, so use the GAP instance inside OSCAR:

```julia
using ArithmeticVoronoi, Oscar
import Oscar: GAP

K = hap_voronoi_complex(voronoi_request(-7; group=:SL);
    resolution_length=2, progress=true);
resolution = GAP.Globals.FreeGResolution(K, 2);
chain = GAP.Globals.TensorWithIntegers(resolution);
sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(chain, 1))))  # [0, 4]: Z × C4
```

This constructs a new complex in Julia and needs no Magma or stored data.
`voronoi_request(D)` takes a field **discriminant**; the older
`hap_voronoi_complex(D)` shortcut still works. HAP's quadratic convenience
constructors take **squarefree radicands**:
their Gaussian case uses `-1`, where our Julia interface uses `-4`.
HAP is not a required component of the native presentation algorithm.

The [HAP and Schönnenbeck guide](docs/hap_and_vmh.md) explains the package
layers, native higher homology, an optional Julia provider for HAP's own
Voronoi entry points, and independent comparisons with fresh Magma complexes.
Native complexes also accept higher rank through `voronoi_request(D; rank=n)`:
quadratic GLₙ and imaginary-quadratic SLₙ on Oⁿ. The presentation pipeline
above remains rank two. See [higher-rank computations and validation](docs/higher_rank.md).
Cup products require additional contracting-homotopy data.

Current comparisons are:

| Computation | Comparison and scope |
| --- | --- |
| 16 real fields | All 898 perfect-form orbits and 5,667 labelled neighbours matched to C++ by exact isometries |
| Five imaginary fields | Perfect-form counts and polytope types agree with Yasaki |
| Five imaginary SL/PSL pairs | Abelianizations agree with HAP's stored complexes |
| Five imaginary GL/SL/PSL triples | Abelianizations agree with Page's Magma code; all exported matrix relations checked in OSCAR |
| D = 5 and −7 | Fresh Schönnenbeck complexes imported into HAP; GL first homology agrees with our abelianizations |
| Native complexes, D = 5 and −7 | Cell counts and GL H₂ agree with the independent VMH results; H₃ also checked for −7 |

For the real and imaginary results, complete tables and timing methodology
are in [the real survey](docs/native_survey.md) and
[the imaginary survey](docs/imaginary_quadratic.md). Matching abelianizations
and matrix relations alone does not certify an isomorphism of presentations.

For D = 5, 8, 12 or 13, `compare_kw(R)` constructs explicit maps for the
Kirchheimer–Wolfart comparison and saves them under `R.directory/KW/`.
To check local copies of the recorded source transcriptions as well, use
`compare_kw(R; reference_source_directory="/path/to/transcriptions")` with
`KWeng.tex` and `KWger.tex`. See [KW_COMPARISON.md](KW_COMPARISON.md) for the
scope, source hashes and D8 transcription discrepancy.

## Optional C++ backend

This section is for C++ comparisons and older scripts. The build, exporter
and pinned upstream sources belong to
[PerfectForms.jl](https://github.com/dimpase/PerfectForms.jl).
The native workflows above do not use them.

### Build and install

The C++ path needs Linux or macOS, Git, CMake 3.24+, Make or Ninja,
pkg-config, a C++20 compiler and standard library with `std::format`, and
development files for GMP, Eigen 3, Boost.Serialization and nauty/Traces.
On macOS install current Xcode Command Line Tools, then:

```sh
brew install cmake ninja pkgconf gmp eigen@3 boost nauty
export CMAKE_PREFIX_PATH="$(brew --prefix eigen@3):$(brew --prefix)"
```

After instantiating this Julia project, build its pinned PerfectForms source
from the checkout root:

```sh
perfect_forms_source="$(julia --startup-file=no --project=. -e 'using PerfectForms; print(pkgdir(PerfectForms))')"
cmake -S "$perfect_forms_source" -B build/perfect-forms \
  -DCMAKE_BUILD_TYPE=Release '-DCMAKE_CXX_FLAGS_RELEASE=-O1 -DNDEBUG' \
  -DCMAKE_INSTALL_PREFIX="$PWD/backend-install"
cmake --build build/perfect-forms --parallel 1
cmake --install build/perfect-forms --prefix "$PWD/backend-install"
export PERFECT_FORMS_BACKEND_CONFIG="$PWD/backend-install/config.json"
```

This also works when Julia's package source is read-only: the build and
installation directories are in this checkout. CMake fetches and verifies
the pinned upstream sources and recursive submodules. Initial configuration
and compilation can take several minutes. The serial `-O1` build reduces
compiler memory use.

For libraries in another prefix, set `CMAKE_PREFIX_PATH`. The nauty overrides
are `Nauty_INCLUDE_DIR` for `traces.h`, `Nauty_CORE_INCLUDE_DIR` for `nauty.h`
and `Nauty_LIBRARY`. System nauty/Traces is required.

### Select the installation

The configuration variable must contain the **absolute path to the installed
`config.json` file**, not its directory. An existing installation can be reused:

```sh
export PERFECT_FORMS_BACKEND_CONFIG="$HOME/work/tmp/PerfectForms.jl/backend/config.json"
```

If Julia is already running, set it inside that session:

```julia
ENV["PERFECT_FORMS_BACKEND_CONFIG"] =
    expanduser("~/work/tmp/PerfectForms.jl/backend/config.json")
Rcpp = run_field_pipeline(12; backend=:polyhedral_common, progress=true,
    directory="results/my_cpp_D12", reuse_baseline=false);
```

One exporter handles all the real survey fields. `config.json` records its
absolute executable path and hash; `backend-build.json` records provenance.
Reinstall after rebuilding or moving it. Changed executable hashes invalidate
C++ enumeration caches. `BackendConfig(absolute_executable_path)` remains
available for explicit configurations.

The older `scripts/run_fields.jl 12` and `scripts/word_example.jl 12` select
the default C++ backend and require this configuration. The native batch
scripts are the recommended commands for Julia-only runs.

The C++ process runner has deadline/memory controls; Linux uses `prlimit`
from util-linux and macOS uses a `ps`-based resident-memory watchdog.
Those controls and sanitizer builds are documented in PerfectForms. They
are separate from the native enumerator. For retained sanitizer exports,
`scripts/validate_sanitizer.jl` (optionally `--fields`) performs our independent
exact checks. See [DEPENDENCIES.txt](DEPENDENCIES.txt) for source attribution.

## Tests

These checks need neither a C++ executable nor the optional archives:

```sh
julia --startup-file=no --project=. test/native.jl
julia --startup-file=no --project=. test/quadratic_units.jl
julia --startup-file=no --project=. test/imaginary.jl
julia --startup-file=no --project=. test/imaginary_bundle.jl
julia --startup-file=no --project=. test/page_comparison.jl
```

The last command replays the committed Page export in OSCAR; it does not
launch Magma. To regenerate that export, follow the imaginary-field guide.
With C++ configured, compare fresh enumerations:

```sh
julia --startup-file=no --project=. test/native_comparison.jl
# Or select a subset:
julia --startup-file=no --project=. test/native_comparison.jl 13 40
```

The full historical regression suite also needs C++ and the D5/D8 fixtures:

```sh
git lfs install --local
git lfs pull
python3 scripts/restore_fixtures.py
julia --startup-file=no --project=. test/runtests.jl
ctest --test-dir build/perfect-forms --output-on-failure
```

Fixture restoration requires Python 3.12+, checks hashes and refuses to
overwrite different data. Tests create local scratch results and can update
tracked reports. The [CI workflow](.github/workflows/ci.yml) builds and tests
the C++ and Julia paths on Linux and macOS.

## Optional historical archives

Git LFS stores the sealed five-field dataset and the original D5/D8 regression
fixtures. They describe the historical C++ delivery, not the current native
survey bundles. Install Git LFS and run `git lfs pull` when you need them.
A small file beginning with `version https://git-lfs` is a pointer whose
archive payload has not yet been downloaded.

[The dataset report](results/FIELDS_WORDS_DATA_REPORT.txt) records the original
counts and timings. Verify and unpack its archive from the checkout root:

```sh
python3 -c 'import hashlib,json,pathlib; p=pathlib.Path("results/dataset"); m=json.loads((p/"archive_manifest.json").read_text()); assert hashlib.sha256((p/m["archive"]).read_bytes()).hexdigest()==m["sha256"]; print("archive checksum OK")'
tar -xzf results/dataset/arithmetic_groups_dataset_v1.tar.gz -C results/dataset
```

The extracted `results/dataset/arithmetic_groups_dataset_v1/` contains field
presentations, maps, geometric and word certificates, `dataset.json`, and the
software snapshot that verified that delivery. Use its own environment:

```sh
julia --startup-file=no --project=results/dataset/arithmetic_groups_dataset_v1/software \
  -e 'using Pkg; Pkg.instantiate(;allow_autoprecomp=false)'
julia --startup-file=no --project=results/dataset/arithmetic_groups_dataset_v1/software \
  results/dataset/arithmetic_groups_dataset_v1/software/scripts/verify_dataset.jl \
  results/dataset/arithmetic_groups_dataset_v1
```

Add `--complete` to replay the geometric construction and Tietze maps.
The default algebraic check verifies matrices, relations, maps, invariants
and supplied words. Neither verification needs a C++ executable. Keep the
extracted bundle unchanged so its manifest remains valid. See
[DATASET_SCHEMA.md](DATASET_SCHEMA.md) for verification levels. The old
`export_dataset` format requires a C++ checkpoint and does not yet accept
native enumerations; use the native survey bundles for those results.

## Troubleshooting and timing

| Symptom | Action |
| --- | --- |
| Long pause at the first `using` or computation | Allow for Julia/OSCAR import and compilation; subsequent calls in the same session are faster |
| `missing_executable` on a native computation | Pass `backend=:native` explicitly; the API default still selects C++ |
| `config.json needs an absolute path` | Set the C++ configuration to an absolute file path, using `expanduser` in Julia if needed |
| Configuration points to `.../backend` | Append `/config.json`; the directory itself is not a configuration file |
| A cached result no longer matches | Use a new output directory; regenerate configuration after changing the C++ executable |
| HAP fails to load | Follow the GAP.jl package installation example in the HAP guide |

Timings in `run.json` exclude process startup but may include compilation
and cache replay. Native enumeration benchmarks in the surveys use warm-up
and fresh enumerations. Cached pipeline runs, enumeration-only timings and
Page's full domain/presentation timings measure different work; compare the
documented scopes before drawing speed conclusions.

Repository code is in `src/`, `scripts/` and `test/`; current run directories
are ignored under `results/`. Historical receipts retain their original paths
and environment. Use this guide for current commands and the
[README](README.md) for a short introduction.
