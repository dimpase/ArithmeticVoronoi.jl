# Validate Page's exported presentations and compare their abelianizations
# with the independently computed native presentations. No Magma is needed
# to replay this check; see docs/imaginary_quadratic.md to regenerate its input.
using ArithmeticVoronoi, Oscar, JSON3, Test
const PCAV = ArithmeticVoronoi
length(ARGS)==2 || error("Usage: julia --project=. test/page_comparison.jl /path/to/page_presentations.json /path/to/imaginary_presentations.json")
page_path, native_path = abspath.(ARGS)
page = JSON3.read(read(page_path, String), Dict{String,Any})
native = JSON3.read(read(native_path, String), Dict{String,Any})

@testset "Page presentations and native abelianizations" begin
    @test page["schema"] == "av_page_presentations_v1"
    expected_keys = Set((D, name) for D in native["discriminants"] for name in ("GL", "SL", "PSL"))
    actual_keys = [(entry["discriminant"], entry["group"]) for entry in page["groups"]]
    @test length(actual_keys) == length(expected_keys)
    @test Set(actual_keys) == expected_keys
    for entry in page["groups"]
        D, name = entry["discriminant"], entry["group"]
        P = arithmetic_problem(D)
        I = identity_matrix(P.field, 2)
        matrices = [PCAV.field_matrix_wire(P, m) for m in entry["generator_matrices"]]
        words = [Int.(w) for w in entry["relators"]]
        @test all(check_group_matrix(P, encode_matrix(P, g)) for g in matrices)
        @test name == "GL" || all(det(g) == 1 for g in matrices)
        @test all(begin
            value = PCAV.evaluate_word(matrices, word)
            value == I || (name == "PSL" && value == -I)
        end for word in words)
        F = free_group(length(matrices))
        G, _ = quo(F, [F(w) for w in words])
        A, _ = maximal_abelian_quotient(G)
        factors = BigInt.(elementary_divisors(abelian_invariants(A)))
        @test sort(factors) == sort(BigInt.(entry["abelian_invariants"]))
        expected = native["fields"][string(D)]["groups"][name]["abelianization"]
        @test count(iszero, factors) == expected["free_rank"]
        @test string.(filter(x -> x > 1, factors)) == expected["torsion"]
    end
end
