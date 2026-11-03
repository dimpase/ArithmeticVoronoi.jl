# Resolve from a fresh client environment, without instantiating the package
# checkout first. --test also exercises Pkg.test and its subprocess environment.
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

if "--test" in ARGS
    Pkg.test("ArithmeticVoronoi";coverage=true)
end
