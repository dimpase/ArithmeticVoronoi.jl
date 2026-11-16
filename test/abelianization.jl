using ArithmeticVoronoi, Oscar, Test, Random

@testset "sparse exact abelianization" begin
    rng=MersenneTwister(8732)
    for m in (0,1,3,8,20),n in (0,1,4,9)
        A=rand(rng,-3:3,m,n)
        words=[reduce(vcat,[fill(sign(A[i,j])*j,abs(A[i,j])) for j in 1:n];init=Int[]) for i in 1:m]
        actual=abelianization_invariants(n,words)
        S=snf(matrix(ZZ,A));diagonal=[abs(S[i,i]) for i in 1:min(m,n) if S[i,i]!=0]
        @test actual.free_rank==n-length(diagonal)
        @test actual.torsion==string.([x for x in diagonal if x>1])
        duplicate=abelianization_invariants(n,[words;[reverse(-w) for w in words]])
        @test (duplicate.free_rank,duplicate.torsion)==(actual.free_rank,actual.torsion)
    end
    n=500;words=[[i,-i-1] for i in 1:n-1];push!(words,[n])
    result=abelianization_invariants(n,repeat(words,20))
    @test result.free_rank==0 && isempty(result.torsion)
    torsion=abelianization_invariants(3,[[1,1,2],[2,2,2]])
    @test torsion.free_rank==1 && torsion.torsion==["6"]
    @test abelianization_invariants(2,[[1,2,-1,-2]]).free_rank==2
    @test_throws ArgumentError abelianization_invariants(2,[[3]])
    @test_throws ArgumentError abelianization_invariants(2,[[0]])
    @test_throws ArgumentError abelianization_invariants(-1,Vector{Int}[])
end
