# Releasing ArithmeticVoronoi

The package is registered from this repository, with UUID
`ad54a9ab-83c5-468c-8b76-196f5a70a2f9`. The first version is `0.1.0`.

1. Update `Project.toml` and check the compatibility bounds. All dependencies
   of the registered engine must be available in General. Keep the optional
   PerfectForms comparison environment separate; it is not a dependency of
   the registered package.
2. Run `julia --startup-file=no test/registry_install.jl --test`. This resolves
   from a fresh client environment, runs a native computation, and runs the
   full suite with coverage without first instantiating the checkout. Require
   the LTS/stable Linux, macOS, HAP and documentation CI jobs to pass on the
   release commit. See [development instructions](../docs/src/development.md)
   for building the manual locally.
3. Enable the [JuliaRegistrator GitHub app](https://github.com/apps/juliateam-registrator)
   for this repository. Comment `@JuliaRegistrator register` on the release
   commit. For a release branch, comment on its commit or use
   `@JuliaRegistrator register branch=BRANCH` in an issue comment.
4. Review the generated General PR and its checks. Correct the package and
   retrigger Registrator if needed, keeping the same version while that
   registration is pending so the existing PR is updated. A new package normally has a three-day
   review period. Registration is complete only after the General PR merges.
5. TagBot creates the tag and GitHub release after registration. No tag or
   GitHub release is needed beforehand. Keep workflow changes out of the final
   release commit; GitHub restricts TagBot's token when tagging such commits.

After registration, test in a new Julia environment:

```julia
using Pkg
Pkg.activate(temp=true)
Pkg.add("ArithmeticVoronoi")
using ArithmeticVoronoi
E = enumerate_perfect_forms(classical_problem())
@assert E.completion_status == :complete
```

ArithmeticVoronoiData has its own repository, UUID and release cycle. Its
registration is independent and does not add a dependency to this engine.

See the [General guidelines](https://github.com/JuliaRegistries/General#registering-a-package-in-general)
and [Registrator instructions](https://github.com/JuliaRegistries/Registrator.jl#how-to-use).
