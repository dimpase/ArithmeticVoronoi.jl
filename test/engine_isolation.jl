using Test

@testset "engine works without reference data or result libraries" begin
    root=normpath(joinpath(@__DIR__,".."))
    mktempdir() do isolated
        # Deliberately copy only runtime source, not test fixtures, references,
        # validation tools, recorded results, or the separate data package.
        cp(joinpath(root,"src"),joinpath(isolated,"src"))
        probe=joinpath(@__DIR__,"engine_isolation_probe.jl")
        # Pkg.test runs in a resolved temporary environment. Reuse it rather
        # than assuming that the package checkout has its own manifest.
        project=dirname(Base.active_project())
        command=`$(Base.julia_cmd()) --startup-file=no --project=$project $probe $isolated`
        @test success(pipeline(command;stdout,stderr))
    end
end
