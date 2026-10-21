using ArithmeticVoronoi, Oscar, JSON3, Test
const IBAV=ArithmeticVoronoi
path=isempty(ARGS) ? joinpath(@__DIR__,"../results/imaginary_presentations.json") : only(ARGS)
bundle=JSON3.read(read(path,String),Dict{String,Any})

@testset "portable imaginary presentations in OSCAR" begin
    @test bundle["schema"]=="av_imaginary_presentations_v1"
    @test bundle["status"]=="complete"
    for D in bundle["discriminants"]
        field=bundle["fields"][string(D)];P=arithmetic_problem(D);I=identity_matrix(P.field,2)
        @test field["comparison"]["status"]=="passed"
        for (name,data) in field["groups"]
            matrices=[IBAV.field_matrix_wire(P,g) for g in data["generator_matrices"]]
            relations=[Int.(w) for w in data["relators"]]
            @test all(check_group_matrix(P,encode_matrix(P,g)) for g in matrices)
            @test name=="GL" || all(det(g)==1 for g in matrices)
            @test all(begin
                value=IBAV.evaluate_word(matrices,w)
                value==I || (name=="PSL" && value==-I)
            end for w in relations)
            F=free_group(length(matrices));G,toG=quo(F,[F(w) for w in relations])
            A,toA=maximal_abelian_quotient(G)
            invariants=BigInt.(elementary_divisors(abelian_invariants(A)))
            expected=data["abelianization"]
            @test count(iszero,invariants)==expected["free_rank"]
            @test string.(filter(x->x>1,invariants))==expected["torsion"]
        end
    end
end
