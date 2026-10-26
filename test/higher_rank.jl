using ArithmeticVoronoi, Oscar, Test
import Oscar: GAP
const HRAV=ArithmeticVoronoi

@testset "higher-rank trace spaces and coordinates" begin
    for n in (2,3,4),D in (1,5,8,-3,-4)
        P=D==1 ? classical_problem(;rank=n) : arithmetic_problem(D;rank=n)
        dimension=D==1 ? div(n*(n+1),2) : (D<0 ? n*n : n*(n+1))
        @test P.rank==n
        @test length(P.form_basis)==dimension
        coefficients=QQ.([i//7 for i in 1:dimension])
        A=HRAV.native_direction(P,coefficients)
        F=form_from_trace(P,A)
        @test trace_gram(P,F)==A
        @test HRAV.form_coordinates(P,F)==coefficients
        g=identity_matrix(P.field,n)
        g[1,n]=P.degree==1 ? one(P.field) : P.generator
        M=encode_matrix(P,g)
        @test decode_matrix(P,M)==g
        @test trace_gram(P,HRAV.act_form(P,F,g))==transpose(M)*A*M
        # The coefficient permutation ceases to be an involution at rank 3.
        @test internal_gram(P,backend_gram(P,A))==A
        z=collect(1:n*P.degree)
        @test internal_vector(P,z[HRAV.backend_permutation(P)])==z
        @test coordinate_vector(P,field_vector(P,z))==QQ.(z)
    end
    @test arithmetic_problem(voronoi_request(-4;rank=3)).rank==3
    @test_throws ArgumentError voronoi_request(-4;rank=1)
    @test_throws ArgumentError classical_problem(;rank=true)
    @test enumerate_perfect_forms(arithmetic_problem(-4;rank=3)).completion_status==:complete
end

@testset "odd-rank special-linear stabilizers retain their actual centre" begin
    P=arithmetic_problem(-3;rank=3)
    I=identity_matrix(P.field,3)
    V=shortest_data(trace_gram(P,I)).vectors
    # Monomial matrices with sixth-root entries: 3!*6^3/6 = 216 in SL.
    stabilizer=HRAV.cell_maps(P,V,V;determinant_one=true)
    @test length(stabilizer)==216
    @test all(g->det(g)==1,stabilizer)
    @test all(P.generator^k*I in stabilizer for k in (0,2,4))
    @test !(-I in stabilizer)
    generators=HRAV.cell_generators(P,stabilizer)
    @test all(g->det(g)==1,generators)
    @test length(HRAV.finite_matrix_closure(encode_matrix.(Ref(P),generators),6,216).elements)==216
end

@testset "classical higher-rank enumeration and homology" begin
    E=enumerate_perfect_forms(classical_problem(;rank=3);backend=:native)
    @test length(E.representatives)==1
    C=voronoi_complex(E)
    @test length.(C.cells)==[1,1,2,1]
    @test validate_voronoi_complex(C)
    R=GAP.Globals.FreeGResolution(hap_complex(C;resolution_length=3),3)
    chain=GAP.Globals.TensorWithIntegers(R)
    @test sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(chain,1))))==[2]
    @test sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(chain,2))))==[2,2]
    S=voronoi_complex(E;group=:SL)
    @test validate_voronoi_complex(S)
    @test all(det(g)==1 for layer in S.stabilizers for stabilizer in layer for g in stabilizer)
    Rsl=GAP.Globals.FreeGResolution(hap_complex(S;resolution_length=2),2)
    @test isempty(GAP.gap_to_julia(GAP.Globals.Homology(GAP.Globals.TensorWithIntegers(Rsl),1)))
    # A4 and D4: independent classical perfect-form and automorphism counts.
    E4=enumerate_perfect_forms(classical_problem(;rank=4);backend=:native)
    @test length(E4.representatives)==2
    @test sort([r.stabilizer_order for r in E4.representatives])==[240,1152]
end

@testset "rank-three Gaussian GL and SL through HAP" begin
    enable_hap_julia_backend(;resolution_length=2)
    construct=GAP.evalstr("r -> VoronoiComplex(r : backend := \"julia\")")
    P=arithmetic_problem(-4;rank=3)
    for (group,expected) in ((:GL,[4]),(:SL,Int[]))
        K=construct(GAP.GapObj(voronoi_request(-4;group,rank=3);recursive=true))
        @test GAP.Globals.IsHapNonFreeResolution(K)
        @test GAP.evalstr("K -> K!.arithmeticVoronoi.rank")(K)==3
        cells=Int.(GAP.gap_to_julia(GAP.evalstr("K -> List([0..6], K!.dimension)")(K)))
        @test cells==[1,1,3,5,4,3,2]
        R=GAP.Globals.FreeGResolution(K,2)
        h1=sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(GAP.Globals.TensorWithIntegers(R),1))))
        @test h1==expected
        if group==:GL
            reference=GAP.Globals.ResolutionGL3QuadraticIntegers(-1,2)
            @test h1==sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(
                GAP.Globals.TensorWithIntegers(reference),1))))
        else
            generators=GAP.gap_to_julia(GAP.evalstr("K -> GeneratorsOfGroup(K!.group)")(K);recursive=true)
            for rows in generators
                M=matrix(QQ,6,6,[rows[i][j] for i in 1:6 for j in 1:6])
                @test det(decode_matrix(P,M))==1
            end
            # -I has determinant -1 in odd rank, and must not enter SL's stabilizers.
            @test GAP.evalstr("K -> ForAll([0..6], d -> ForAll([1..K!.dimension(d)], j -> not -IdentityMat(6) in K!.stabilizer(d,j)))")(K)
        end
    end
end
