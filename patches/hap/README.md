# Optional Julia provider for HAP's Voronoi functions

`0001-optional-voronoi-backends.patch` is prepared against
[`gap-packages/hap` at `91bd7baad0ef56c91022cbd6ba2a709942ce2e15`](https://github.com/gap-packages/hap/tree/91bd7baad0ef56c91022cbd6ba2a709942ce2e15).
It loads the Voronoi entry points at package startup, adds a general record
interface and provider registration, and includes a GAP test and tutorial
documentation. Loading these definitions does not run Magma. It is submitted
upstream as [draft HAP PR #167](https://github.com/gap-packages/hap/pull/167).

In a HAP checkout, apply it with:

```sh
git apply /absolute/path/to/arithmetic_voronoi/patches/hap/0001-optional-voronoi-backends.patch
```

The interface is independent of Julia. An external package registers a
function `provider(request)` using
`RegisterVoronoiBackend("name", provider)`. HAP passes a request record
unchanged; providers interpret the field, order, lattice, group and any
additional data they support. The provider returns a HAP
non-free resolution with full stabilizers, orientations and boundaries.
The caller chooses it explicitly:

```gap
request := rec(group:="GL", rank:=2,
    field:=rec(kind:="quadratic", discriminant:=-7),
    order:=rec(kind:="maximal"), lattice:=rec(kind:="free"));
K := VoronoiComplex(request : backend := "julia");
R := FreeGResolution(K,4);
List([1..3], n -> Homology(TensorWithIntegers(R),n));
```

The general constructor requires an explicit registered backend. The legacy
`VoronoiComplexGL(n,d)` and `VoronoiComplexSL(n,d)` calls without a backend
option retain the Magma path; with a provider they construct quadratic records
using the squarefree radicand `d`. An unavailable provider
is an error; an unsupported request never silently launches another backend.
HAP acquires no Julia dependency. The patch test uses a tiny mock provider
and requires neither Julia nor Magma:

```gap
LoadPackage("hap");
Test(Filename(DirectoriesPackageLibrary("hap","tst/testextra"),
              "voronoi-backends.tst"));
```

ArithmeticVoronoi registers its real provider inside Julia/GAP.jl:

```julia
using ArithmeticVoronoi, Oscar
import Oscar: GAP
enable_hap_julia_backend()
construct = GAP.evalstr("r -> VoronoiComplex(r : backend := \"julia\")");
K = construct(GAP.GapObj(voronoi_request(-7); recursive=true));
```

For released HAP without this patch, `enable_hap_julia_backend()` loads
[a compatibility wrapper](../../gap/VoronoiBackends.g) implementing the same
dispatch. Installed HAP files are not modified. The callback runs in the
existing Julia process; standalone GAP would need its own registered
provider or a Julia-enabled environment.

See [request records and experiment runners](../../docs/requests.md) for the
shared Julia/JSON interface, and [the HAP guide](../../docs/hap_and_vmh.md) for supported groups and the
direct `hap_voronoi_complex` Julia API. Native integration tests run with:

```sh
julia --startup-file=no --project=. -e 'include("test/hap.jl")'
```

The integration test is separate from the default suite because HAP is an
optional GAP package. It compares fresh native complexes and higher homology
against the independent Magma/VMH results for discriminants 5 and −7, checks
imaginary SL against a stored HAP resolution, and exercises the provider.

Validated on 2026-10-21 with Julia 1.13.2, OSCAR 1.8.2 and GAP 4.16.1:

- All 55 integration assertions passed with installed HAP 1.78 and the
  compatibility wrapper, in a fresh Julia process.
- All 55 passed with the full patched upstream checkout (which declares
  version 1.80), including its package-startup registration. The GAP provider
  dispatch test also passed.
- All 41 request-validation assertions and 1,185 assertions in `test/native.jl`
  and `test/imaginary.jl` passed, including the record-based presentation pipeline.
- GAPDoc parsed and validated the complete updated tutorial, including the
  general request example.
- The presentation runners completed fresh JSON-request batches for D=5 and
  D=-7, including word certificates and the imaginary comparisons. The HAP
  runner computed H1/H2 for GL at D=5,-7 and SL at D=-7; all ten runner and
  resume checks passed.
- The patch applies cleanly to the pinned unmodified upstream tree.
