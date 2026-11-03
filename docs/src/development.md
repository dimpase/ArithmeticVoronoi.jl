# Development and documentation

The package supports Julia 1.10 and later. The repository contains compatibility
bounds rather than a committed root manifest. Applications and reproducible
experiments should keep their own resolved environments.

## Tests and coverage

CI resolves the package in a fresh client environment and calls
`Pkg.test("ArithmeticVoronoi"; coverage=true)` on Julia LTS and the current
stable release. It also tests macOS and bundled/upstream HAP. The native suite
requires neither the data repository nor a C++ build.

To reproduce the client installation test from a checkout:

```sh
julia --startup-file=no test/registry_install.jl --test
```

CI processes Julia coverage files into `lcov.info` and saves the report as an
Actions artifact for each Julia/platform combination. After all native and HAP
checks pass on `main`, CI uploads the Linux stable report to
[Codecov](https://codecov.io/gh/dimpase/ArithmeticVoronoi.jl) using GitHub OIDC.
No stored Codecov token is needed. Pull requests collect artifacts without
uploading to Codecov.

## Build the manual locally

From the repository root:

```sh
julia --startup-file=no --project=docs -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --startup-file=no --project=docs docs/make.jl
```

Open `docs/build/index.html`. `docs/Manifest.toml` and `docs/build/` are ignored.
The build checks exported docstrings and internal links. The documentation
workflow builds pull requests and `main`, and saves the generated HTML as a
GitHub Pages artifact. Successful builds on `main` publish the
[online manual](https://dimpase.github.io/ArithmeticVoronoi.jl/).
Pull request builds do not deploy.

The historical comparisons record exact software versions in their own
validation receipts; those versions are evidence about those runs, not pins
on subsequent package installations.
