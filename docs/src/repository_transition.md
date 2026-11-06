# Repository name transition

The current engine repository is
[dimpase/ArithmeticVoronoi.jl](https://github.com/dimpase/ArithmeticVoronoi.jl).
It was copied from `dimpase/arithmetic_voronoi`, preserving both branches
(`main` and `m0`), all commit IDs, and the historical LFS objects. The copy
starts at engine commit `14d84b03a3cdfc5776d25ffd0283b78ec3df7a7f`; subsequent
documentation changes use the new name. Current `main` has no LFS-managed
files. The independently maintained data package is
[ArithmeticVoronoiData.jl](https://github.com/dimpase/ArithmeticVoronoiData.jl).

The original repository remains available while HAP updates its reference.
HAP's `.github/workflows/CI.yml` currently selects:

```yaml
repository: dimpase/arithmetic_voronoi
ref: 007df75debde04c7469124891ac78d555d5de20b
```

[HAP PR #169](https://github.com/gap-packages/hap/pull/169), submitted after
ArithmeticVoronoi v0.1.0 was registered, updates HAP to
`dimpase/ArithmeticVoronoi.jl` at the published provider revision `c6942b9` and
expands the Julia examples to rational/quadratic GL, SL and PSL.
It also documents the direct `hap_voronoi_complex` calls. The replacement
patch and application instructions are in [patches/hap/README.md](https://github.com/dimpase/ArithmeticVoronoi.jl/blob/hap-interface/patches/hap/README.md).
The earlier URL-only patch is retained as historical context.

The provider revision is public and the HAP PR is open. Before retiring the
old repository, merge the follow-up in HAP and verify its Julia CI job using
the new location. The new rational request and projective-complex API extends
registered v0.1.0, so the HAP tutorial installs the tested source revision
explicitly.

Existing working copies can follow the new repository without changing their
directory name:

```sh
git remote set-url origin git@github.com:dimpase/ArithmeticVoronoi.jl.git
git fetch origin
```

The Julia package remains `ArithmeticVoronoi`, with the same UUID and
`using ArithmeticVoronoi` entry point. The old repository has not been
renamed, deleted or archived, so existing HAP references continue to work.
