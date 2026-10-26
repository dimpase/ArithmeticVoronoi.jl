# Exercise package installation as a dependency: the checkout's Manifest.toml
# and the separate C++ comparison environment must not affect resolution.
using Pkg, Test

Pkg.activate(mktempdir())
Pkg.develop(path=normpath(joinpath(@__DIR__,"..")))
Pkg.instantiate(;allow_autoprecomp=false)

using ArithmeticVoronoi

@testset "installation from registered dependencies" begin
    dependencies=values(Pkg.dependencies())
    @test !any(d->d.name=="PerfectForms",dependencies)
    @test !any(d->d.name=="ArithmeticVoronoiData",dependencies)
    mktempdir() do output
        cd(output) do
            R=run_field_pipeline(-4;simplify=false)
            @test R.run["backend"]=="native"
            @test !R.enumeration.provenance["from_cache"]
            @test realpath(R.directory)==joinpath(realpath(output),"results","fields","D-4")
            @test isfile(joinpath(R.directory,"run.json"))
        end
    end
end
