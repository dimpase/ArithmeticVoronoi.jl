# Arithmetic computation requests

Perfect-form enumeration, presentation experiments and the Julia HAP provider
accept the same record describing the group, rank, field, order and lattice:

```julia
using ArithmeticVoronoi

request = (
    group = "GL", rank = 2,
    field = (kind="quadratic", discriminant=-7),
    order = (kind="maximal",),
    lattice = (kind="free",),
);
# Equivalent convenience constructor:
request = voronoi_request(-7);

E = enumerate_perfect_forms(request; progress=true);
R = run_field_pipeline(request; progress=true,
    directory="results/request_D-7", reuse_baseline=false);
K = hap_voronoi_complex(request; resolution_length=3);
```

These record-based calls use native Julia enumeration by default. The legacy
`run_field_pipeline(D)` and problem-object enumeration calls retain their
existing default; pass `backend=:native` to them explicitly. The presentation
runner `run_field_pipeline` takes a rank-two GL request and computes GL, SL
and PSL presentations; `run_gl_pipeline` accepts any supported rank and
computes GL. HAP
complexes currently support GL and imaginary-quadratic SL. For the latter,
use `voronoi_request(-7; group=:SL)`.

Julia accepts named tuples, dictionaries with string or symbol keys, and
JSON records. `normalize_voronoi_request(request)` validates the supported
arithmetic and returns a dictionary suitable for JSON and GAP.jl. A quadratic
field can specify either a fundamental `discriminant` or a squarefree
`radicand`, never both. For example, `radicand=-1` normalizes to
`discriminant=-4`.

Native enumeration, `run_gl_pipeline` and HAP complexes accept any integer rank at least two,
a quadratic field, its maximal order, and the free lattice `O^rank`.
For example, `voronoi_request(-4; rank=3)` requests GL₃ over the Gaussian
integers. Other fields, nonmaximal orders, nonfree lattices and unknown record
fields are rejected before enumeration. The lower-level
`classical_problem(; rank=3)` API supports rational forms as well; the request
adapter covers quadratic fields. The C++ adapter and the GL/SL/PSL presentation
pipeline still require rank two. See [higher-rank validation](higher_rank.md).

```julia
R = run_gl_pipeline(voronoi_request(-4; rank=3);
    directory="results/request_GL3_D-4");
```

The GL-only batch runner uses the same records:

```sh
julia --project=. scripts/run_gl.jl --requests=examples/higher_rank_gl_requests.json --output=results/higher_gl
```

## Run experiments from JSON

[examples/voronoi_requests.json](../examples/voronoi_requests.json) contains
GL requests for discriminants -7 and 5. From the checkout root:

```sh
# GL/SL/PSL presentations and matrix-to-word certificates:
julia --project=. scripts/run_native_fields.jl --requests=examples/voronoi_requests.json --output=results/request_presentations

# New complexes and integral homology H_1 and H_2 through GAP/HAP:
julia --project=. scripts/run_hap.jl --requests=examples/voronoi_requests.json --degree=2 --output=results/request_homology
```

The imaginary presentation runner also accepts `--requests=FILE`; supply only
negative discriminants and GL requests. It retains its independent comparisons
against HAP's stored SL/PSL complexes. `run_hap.jl` instead constructs fresh
complexes using the Julia provider and accepts imaginary SL requests too.
Its `summary.json` and individual group files contain the request, timings,
and homology invariant factors (`0` denotes an infinite cyclic summand).
The expanded [rank-two HAP batch](../examples/hap_requests.json) adds SL,
Gaussian units, and a GL example at discriminant −15; pass it with
`--requests=examples/hap_requests.json` to the same runner.

For rank-three Gaussian GL and SL in one Julia session:

```sh
julia --project=. scripts/run_hap.jl --requests=examples/higher_rank_hap_requests.json --degree=1 --output=results/rank3_homology
```

Higher-rank filenames include the rank, for example `GL3_D-4.json`, so they
cannot overwrite a rank-two result in the same batch.

Integer arguments remain shorthand for quadratic GL requests, so existing
batch commands still work. Both types of runner accept `--resume`; a saved
result is retained only if its complete request matches. The homology runner
also checks the requested degree. Older results without a request are
recomputed. Presentation runs keep their existing history and certificates.

## General HAP interface

Inside Julia, register the provider and call HAP with the record:

```julia
using ArithmeticVoronoi, Oscar
import Oscar: GAP
enable_hap_julia_backend(; resolution_length=3)
construct = GAP.evalstr("r -> VoronoiComplex(r : backend := \"julia\")");
K = construct(GAP.GapObj(voronoi_request(-7); recursive=true));
R = GAP.Globals.FreeGResolution(K, 3);
GAP.gap_to_julia(GAP.Globals.Homology(
    GAP.Globals.TensorWithIntegers(R), 2)) # [2,12]
```

The interface merged in [HAP PR #167](https://github.com/gap-packages/hap/pull/167)
registers callbacks of the form `provider(request)`. HAP passes the record
unchanged, allowing other providers to interpret field polynomials, order
bases, lattice ideals or subgroup data. Validation of those data belongs to
the provider. The two existing quadratic GL/SL constructors are convenience
wrappers when a provider is selected; their arguments remain radicands.

`enable_hap_julia_backend()` preserves upstream HAP's entry points and provides
the same interface on older HAP releases via a compatibility wrapper.
All computation stays in the Julia/GAP process;
no Magma process or intermediate export file is used. HAP still requires the
stabilizers, orientation actions and equivariant boundaries described in
[the HAP guide](hap_and_vmh.md). The record interface does not add a
contracting homotopy or change which cohomology operations are available.
