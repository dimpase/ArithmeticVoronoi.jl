# ArithmeticVoronoi

Compute presentations and verified matrix-to-word expressions for **GL(n,O)**
(`n ≥ 2`), and **SL(2,O), PSL(2,O)**, where O is the maximal order of a real or
imaginary quadratic field, acting on the free lattice Oⁿ. GL(n,Z) is also
supported. The native backend uses
Julia, OSCAR and Hecke; no C++ build or Magma installation is needed.

The implemented group algorithm is from Oliver Braun, Renaud Coulangeon,
Gabriele Nebe and Sebastian Schönnenbeck,
[*Computing in arithmetic groups with Voronoi's algorithm*](https://arxiv.org/abs/1407.6234):
perfect forms, their neighbours and stabilizers, and ridge relations give a
finite presentation and constructive membership algorithm.

## Start a computation

Use Julia 1.13 (tested with 1.13.2). From a new checkout:

```sh
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/dimpase/arithmetic_voronoi.git
cd arithmetic_voronoi
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(;allow_autoprecomp=false)'
julia --startup-file=no --project=.
```

At the Julia prompt, compute all three groups for O = Z[√3] and work with SL:

```julia
using ArithmeticVoronoi, Oscar

R = run_field_pipeline(12; backend=:native, progress=true,
    directory="results/my_D12", reuse_baseline=false);
S = R.simplified[:SL];
F = free_group(length(S.generator_matrices));
G, toG = quo(F, [F(r) for r in S.relators]);
A, toA = maximal_abelian_quotient(G);
elementary_divisors(abelian_invariants(A))   # [2, 6]: C2 × C6
```

The input is a **field discriminant**: `12` means Q(√3), `-7` means Q(√−7),
and `-4` means Q(i). Change the discriminant and output directory together.
Select `:GL` or `:PSL` to use the other computed presentations. GL and SL
retain −I; PSL is SL/{I,−I}.

Pass `backend=:native` explicitly: the API default still selects the optional
C++ backend. The first Julia import and compilation can take several minutes.
Saved archives and Python are optional.

## Documentation and results

- [User guide](USER_GUIDE.md): installation, group operations, matrix words,
  batches, saved results, optional C++ backend and troubleshooting.
- [Native algorithm](docs/native_backend.md) and
  [presentation justification](PRESENTATION_JUSTIFICATION.md).
- [16 real quadratic fields](docs/native_survey.md): full presentations and
  exact comparisons with the C++ enumerator.
- [Five imaginary quadratic fields](docs/imaginary_quadratic.md): full
  presentations and comparisons with Yasaki, HAP and Aurel Page's Magma code.
- [Native higher homology with HAP and Schönnenbeck comparisons](docs/hap_and_vmh.md).
- [Higher-rank GL presentations and HAP complexes](docs/higher_rank.md), with Magma regression comparisons.
- [Shared request records and JSON experiment batches](docs/requests.md).

See [LICENSE](LICENSE) and [dependency attribution](DEPENDENCIES.txt).
