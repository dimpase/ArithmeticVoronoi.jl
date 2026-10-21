# HAP from Julia and Schönnenbeck's Voronoi software

## How the packages fit together

HAP (Homological Algebra Programming) is a **GAP package**. You can call it
from Julia using **GAP.jl**, which runs GAP within the Julia process. OSCAR
already includes GAP.jl and its GAP runtime. Loading HAP does not launch
Magma or build a C++ backend.

| Layer | Role | How it is selected |
| --- | --- | --- |
| ArithmeticVoronoi | Native perfect forms, presentations, matrix words and full cell complexes | `using ArithmeticVoronoi` |
| OSCAR | Julia interface to algebra systems, including GAP and Hecke | `using Oscar` |
| GAP.jl | Calls GAP functions and handles GAP objects in Julia | `import Oscar: GAP` |
| HAP | Resolutions, homology and related constructions in GAP | `GAP.Packages.load("hap")` |
| Magma + VMH | Optional independent construction of new arithmetic cell complexes | Run Magma separately, then read its output through GAP.jl |

HAP is managed by GAP's package system, rather than by `Pkg.add("HAP")` or
`using HAP`. Outside this project you can install the Julia package `GAP`
with Julia's package manager and use `using GAP` directly. See the
[GAP.jl introduction](https://github.com/oscar-system/GAP.jl) and
[package-management documentation](https://oscar-system.github.io/GAP.jl/stable/packages/).

## Compute a new complex in Julia, then its homology in HAP

From the repository root, start `julia --startup-file=no --project=.`:

```julia
using ArithmeticVoronoi, Oscar
import Oscar: GAP

# Fresh GL(2,O) computation for Q(sqrt(-7)); input is a field discriminant.
K = hap_voronoi_complex(-7; resolution_length=4, progress=true);
R = GAP.Globals.FreeGResolution(K, 4);
C = GAP.Globals.TensorWithIntegers(R);
[sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(C,n)))) for n in 1:3]
# [[2,2,2], [2,12], [2,2,2,2,12]]
```

This starts with native perfect-form enumeration. It constructs all bounded
Ryshkov cells, finite stabilizers, their orientation actions, and signed
boundary incidences. It uses neither Magma nor HAP's stored arithmetic
complexes. Julia checks the equivariant boundary equation `d²=0` exactly.
HAP then resolves the finite cell stabilizers to the requested degree.
`resolution_length=n+1` suffices for computing Hₙ.

Supported groups are **GL₂ over real and imaginary quadratic maximal orders**
and **SL₂ over imaginary quadratic maximal orders**, on the free lattice O².
Use `group=:SL` for the latter. Central −I is retained. This interface does
not currently supply real-quadratic SL, projective quotients, rank 3 or
nonfree lattices.

To reuse a perfect-form enumeration or inspect the native cells:

```julia
E = enumerate_perfect_forms(arithmetic_problem(-7); backend=:native);
V = voronoi_complex(E; group=:GL);
length.(V.cells)                  # [1,2,1], starting in dimension zero
validate_voronoi_complex(V)       # true
K = hap_complex(V; resolution_length=4);
```

Each HAP object owns its data, so constructing another complex does not
overwrite the first. The returned non-free resolution and the standard
`FreeGResolution` path have no contracting homotopy. Additive homology and
cohomology work; HAP's integral cup-product routine needs additional data.

### Use HAP's Voronoi entry points with the Julia provider

```julia
enable_hap_julia_backend(; progress=true)
K = GAP.evalstr("VoronoiComplexGL(2,-7 : backend := \"julia\")");
S = GAP.evalstr("VoronoiComplexSL(2,-7 : backend := \"julia\")");
```

These HAP calls take **squarefree radicands**: `-1` means Q(i), whose
discriminant in `hap_voronoi_complex` is `-4`. Unsupported Julia requests
raise an error; they do not launch Magma as a fallback. The callback runs
inside the existing Julia/GAP.jl process, with no subprocess or intermediate
file. Existing calls without `backend := "julia"` retain HAP's behavior.

The [HAP patch](../patches/hap/README.md) introduces a general optional
provider interface. `enable_hap_julia_backend()` also loads a compatibility
wrapper for unpatched HAP, so no installed package files need editing.

## Use complexes already stored in HAP

From the repository root, start `julia --startup-file=no --project=.` and run:

```julia
using Oscar
import Oscar: GAP

@assert GAP.Packages.load("hap")

# SL(2,O) for Q(sqrt(-7)), using a complex already stored in HAP.
R = GAP.Globals.ResolutionSL2QuadraticIntegers(-7, 2);
C = GAP.Globals.TensorWithIntegers(R);
sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(C, 1))))   # [0, 4]

# The corresponding PSL group.
Rpsl = GAP.Globals.ResolutionPSL2QuadraticIntegers(-7, 2);
sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(
    GAP.Globals.TensorWithIntegers(Rpsl), 1))))            # [0, 2]

# A stored real quadratic GL example: Q(sqrt(5)).
Rgl = GAP.Globals.ResolutionGL2QuadraticIntegers(5, 2);
sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(
    GAP.Globals.TensorWithIntegers(Rgl), 1))))             # [0, 2]
```

These are invariant factors of H₁ with integral coefficients: `[0,4]` means
Z × C4, while `[0,2]` means Z × C2. The resolution length is 2 because H₁
needs the differential from degree 2. These three examples were executed
successfully on this host.

HAP's arguments here are **squarefree radicands**, not our field
discriminants: use `-1` for Q(i) and `3` for Q(√3), corresponding to our
`-4` and `12`. These stored-complex interfaces support a finite collection
of fields; they are not generic new enumerations for arbitrary inputs.
[HAP's tutorial](https://gap-packages.github.io/hap/tutorial/chap14.html)
also describes its separate Bianchi/Swan computations and the provenance
of its stored complexes. The existing imaginary survey calls the
stored-complex SL/PSL interfaces above.

GAP values remain GAP objects until converted. `Int` arguments work directly;
use `GAP.GapObj("text or a path")` when passing Julia strings to a GAP
function, and `GAP.gap_to_julia` to retrieve an ordinary Julia value.

### If HAP is missing, and how to check its version

Try loading it first. If `load` returns `false`, install it with GAP.jl:

```julia
@assert GAP.Packages.install("hap"; interactive=false)
@assert GAP.Packages.load("hap")
```

The installation is managed in GAP.jl's package area; no system GAP
installation is needed. Check the package actually loaded in this session:

```julia
GAP.gap_to_julia(GAP.Globals.GAPInfo.PackagesLoaded.hap[2])
GAP.Packages.locate_package("hap")
```

The checks recorded here used Julia 1.13.2, OSCAR 1.8.2, GAP.jl 0.17.6 and
HAP 1.78. Older survey receipts record HAP 1.70. Julia's manifest identifies
GAP.jl, but GAP also searches its own available packages; a separately
installed GAP package can affect the HAP version selected. Record the loaded
version when comparing computations.

## Schönnenbeck's packages and the HAP connection

Two relevant repositories provide Magma algorithms and GAP readers:

- [VMH-TotallyRealNumberFields](https://github.com/schoennenbeck/VMH-TotallyRealNumberFields)
  handles the real quadratic examples within its totally real field scope.
- [VMH-ImaginaryQuadraticNumberFields](https://github.com/schoennenbeck/VMH-ImaginaryQuadraticNumberFields)
  is the companion for imaginary quadratic fields, developed by Oliver Braun
  and Sebastian Schönnenbeck.

Their common workflow is:

```text
Magma: VoronoiAlgorithm → Presentation / ComputeComplexGL
                                      ↓ exported complex
Julia → GAP.jl → HAP: ResolutionFromWellRoundedComplex
                     → FreeGResolution → TensorWithIntegers → Homology
```

The exported complex includes cell orbits, stabilizers, orientation-preserving
stabilizers and boundary incidences with group elements. The GAP reader
constructs a non-free resolution object; HAP turns this into a free resolution
to the requested degree. This is the connection described in Schönnenbeck's
[*Resolutions for unit groups of orders*](https://arxiv.org/abs/1609.08835).

Both routes were tested with Magma 2.30-1 and HAP inside Julia. Only the
documented `Magma/BasicData.m` input was changed; algorithms and GAP readers
were used unchanged. The [comparison receipt](../results/vmh_comparison.json)
records revisions, complex hashes, versions and results:

| Field discriminant | Fresh complex cell orbits by dimension | GL first homology | Native comparison |
| ---: | --- | --- | --- |
| 5 | 2, 4, 5, 4, 1 | Z × C2 | Matched |
| −7 | 1, 2, 1 | C2 × C2 × C2 | Matched |

This checks the fresh export/import path and H₁. It does not certify an
isomorphism between the presentations or check higher homology.

The subsequent native integration tests in [test/hap.jl](../test/hap.jl)
also reproduce the independent VMH cell counts and GL H₂:
`[2,12]` for −7 and `[2,2,2,2]` for 5. For −7, GL H₃ is
`[2,2,2,2,12]`. Native SL H₁ and H₂ for −7 agree with HAP's stored
resolution. GL and SL H₁ also agree with the presentation results for
−3, −4, −15 and −23. These checks concern the cell complex and its homology;
they do not supply a contracting homotopy or a cohomology ring.

### Reproduce a fresh real quadratic example

From the ArithmeticVoronoi checkout root, clone the tested VMH revision to
a fresh working directory:

```sh
AV_REPO="$PWD"
export VMH_ROOT="$(mktemp -d /tmp/av-vmh.XXXXXX)"
git clone https://github.com/schoennenbeck/VMH-TotallyRealNumberFields.git "$VMH_ROOT"
git -C "$VMH_ROOT" checkout befcfeda041781b7de922a16a37388decee2c642
cat > "$VMH_ROOT/Magma/BasicData.m" <<'MAGMA'
K := QuadraticField(5);
n := 2;
V := KMatrixSpace(K,n,1);
L := [IntegralBasis(K)[i]*V.j : i in [1..Degree(K)], j in [1..n]];
MAGMA
```

With a licensed `magma` command on PATH, construct the presentation and complex:

```sh
cd "$VMH_ROOT"
magma -b <<'MAGMA'
AttachSpec("Magma/VMH-spec");
SetSeed(1);
V := VoronoiAlgorithm();
G, g := Presentation(V);
AbelianQuotient(G);                 // Z/2 + Z for the real D5 example
ComputeComplexGL("complex.g", V);
quit;
MAGMA
cd "$AV_REPO"
julia --startup-file=no --project=.
```

At the Julia prompt, read that newly constructed complex and compute H₁:

```julia
using Oscar
import Oscar: GAP
@assert GAP.Packages.load("hap")

vmh = ENV["VMH_ROOT"];
GAP.Globals.Read(GAP.GapObj(joinpath(vmh, "GAP", "ReadWellRoundedComplex.gi")));
Kcomplex = GAP.Globals.ResolutionFromWellRoundedComplex(
    GAP.GapObj(joinpath(vmh, "complex.g")), 2);
R = GAP.Globals.FreeGResolution(Kcomplex, 2);
sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(
    GAP.Globals.TensorWithIntegers(R), 1))))   # [0, 2] for D5
```

The reader is in the actual repository's `GAP/` directory, not the illustrative
`gapstuff/` path in its upstream README. `BasicData.m` is under `Magma/`.
Use a new Magma process when changing that input: the package imports its
configuration and caches mathematical data.

### Imaginary quadratic variant

For the tested Q(√−7) case, replace the first clone/configuration block with:

```sh
AV_REPO="$PWD"
export VMH_ROOT="$(mktemp -d /tmp/av-vmh-imaginary.XXXXXX)"
git clone https://github.com/schoennenbeck/VMH-ImaginaryQuadraticNumberFields.git "$VMH_ROOT"
git -C "$VMH_ROOT" checkout b5b79a8103c648c36e90b8bd95880782a3ff0d2e
cat > "$VMH_ROOT/Magma/BasicData.m" <<'MAGMA'
n := 2;
d := -7;
steinitz := 1;
MAGMA
```

Run the same Magma and Julia blocks in fresh processes. Both the Magma
abelianization and the HAP H₁ now give `[2,2,2]`. Here `d` is a radicand and
`steinitz=1` selects the free lattice. The package has subgroup options, but
this check used full GL. Its `projective` presentation option concerns the
quotient by −I; do not identify that with quotienting GL by all scalar units.

## Native complex construction

The native construction uses the bounded faces of the relative Ryshkov
polyhedron. A cell is determined by its complete minimum set S; its
dimension is the dimension of the form space minus the rank of the evaluation
vectors for S. Faces whose minimum sets do not span K² are discarded.
All local cone faces are obtained by intersecting facets, then identified
under the arithmetic group.

In rank two, all transporters between spanning minimum sets can be found
exactly by enumerating the images of two K-independent vectors. This also
computes complete stabilizers and every boundary embedding. The orientation
character is the determinant sign on the cell's rational tangent space;
boundary signs use the outward-normal-first convention. The validation checks
incidences, orientations and cancellation of `d²` in the signed permutation
modules, before any tensoring with the integers.

For imaginary SL, the finite determinant quotient splits the GL cell orbits,
and the boundary maps and stabilizers are restricted to determinant one.
This is the actual SL group, including its centre. No presentation-to-resolution
shortcut is assumed.

For H₁ alone, abelianizing a native presentation remains sufficient. The full
complex provides the extra geometric data needed for higher degrees.
