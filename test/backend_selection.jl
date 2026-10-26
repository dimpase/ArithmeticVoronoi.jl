using ArithmeticVoronoi, Test

@testset "native default and explicit optional backend" begin
    P=classical_problem()
    E=enumerate_perfect_forms(P)
    @test E.completion_status==:complete
    @test !E.provenance["from_cache"]
    @test !isdefined(ArithmeticVoronoi,:PerfectForms)
    @test_throws ArgumentError enumerate_perfect_forms(P;backend=:unknown)
    @test_throws r"integrations/perfect_forms/adapter.jl" enumerate_perfect_forms(P;backend=:polyhedral_common)
    mktempdir() do directory
        path=joinpath(directory,"nested","enumeration.json")
        @test write_enumeration(path,E)==path
        @test isfile(path)
        @test length(ArithmeticVoronoi.file_sha(path))==64
    end
end
