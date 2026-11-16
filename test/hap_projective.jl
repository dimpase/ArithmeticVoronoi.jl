using ArithmeticVoronoi, Oscar, Test
import Oscar: GAP

@testset "rational requests and full-centre PSL complexes" begin
    enable_hap_julia_backend(;resolution_length=3)
    # The adapter corrects the old quotient inverse type without editing HAP.
    GAP.Globals.Read(GAP.GapObj(joinpath(ArithmeticVoronoi.PACKAGE_ROOT,"gap","ArithmeticVoronoi.g")))
    @test GAP.evalstr("IsHapQuotientElement(InverseMutable(QuotientGroup((1,2,3,4,5,6),[(),(1,4)(2,5)(3,6)])))")
    @test GAP.evalstr("Size(DirectProduct(Group(QuotientGroup((1,2,3,4,5,6),[(),(1,4)(2,5)(3,6)])),Group(QuotientGroup((1,2,3,4,5,6),[(),(1,4)(2,5)(3,6)]))))") == 9
    construct=GAP.evalstr("r -> VoronoiComplex(r : backend := \"julia\")")
    homology(K,n)=sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(
        GAP.Globals.TensorWithIntegers(GAP.Globals.FreeGResolution(K,n+1)),n))))
    # Known rational homology: PSL2(Z) = C2 * C3, SL3(Z) is perfect.
    for (group,n,expected) in ((:GL,3,[2]),(:SL,3,Int[]),(:PSL,3,Int[]),(:PSL,2,[6]))
        K=construct(GAP.GapObj(voronoi_request(;group,rank=n);recursive=true))
        @test GAP.Globals.IsHapNonFreeResolution(K)
        @test Int(GAP.evalstr("K -> K!.arithmeticVoronoi")(K).rank)==n
        @test GAP.evalstr("K -> K!.arithmeticVoronoi")(K).groupKind==GAP.GapObj(string(group))
        @test Int(GAP.evalstr("K -> K!.arithmeticVoronoi")(K).discriminant)==1
        @test homology(K,1)==expected
        if group==:PSL
            @test Int(GAP.evalstr("K -> K!.arithmeticVoronoi")(K).centerOrder)==gcd(n,2)
        end
        if group==:PSL && n==2
            @test homology(K,2)==Int[]
        end
    end
    # Fresh direct calls need neither a GAP record nor provider registration.
    K=hap_voronoi_complex(;group=:SL,resolution_length=3)
    reference=GAP.Globals.ResolutionSL2Z(1,3)
    @test homology(K,1)==[12]
    @test homology(K,2)==sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(
        GAP.Globals.TensorWithIntegers(reference),2))))
    for (D,d) in ((-3,-3),(-4,-1),(-7,-7))
        S=hap_voronoi_complex(arithmetic_problem(D);group=:SL,resolution_length=3)
        Q=construct(GAP.GapObj(voronoi_request(D;group=:PSL);recursive=true))
        reference=GAP.Globals.ResolutionPSL2QuadraticIntegers(d,3)
        @test Int(GAP.evalstr("K -> K!.arithmeticVoronoi")(Q).centerOrder)==2
        @test GAP.evalstr("K -> K!.arithmeticVoronoi")(Q).groupKind==GAP.GapObj("PSL")
        for n in 1:2
            @test homology(Q,n)==sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(
                GAP.Globals.TensorWithIntegers(reference),n))))
        end
        # Every cell stabilizer loses precisely the central kernel.
        @test GAP.evalstr("{S,Q} -> ForAll([0..2], n -> ForAll([1..S!.dimension(n)], j -> Size(S!.stabilizer(n,j)) = 2*Size(Q!.stabilizer(n,j))))")(S,Q)
    end
    @test_throws ArgumentError hap_voronoi_complex(;group=:PGL)
    @test_throws ArgumentError hap_voronoi_complex(arithmetic_problem(5);group=:PSL)
    @test_throws GAP.GAPError construct(GAP.GapObj(voronoi_request(5;group=:PSL);recursive=true))
end
