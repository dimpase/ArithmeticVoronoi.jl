using ArithmeticVoronoi, Oscar, Test
import Oscar: GAP

@testset "certified real-quadratic SL2 provider" begin
    enable_hap_julia_backend(;resolution_length=4)
    selected=get(ENV,"ARITHMETIC_VORONOI_HAP_ROOT","")
    isempty(selected) || @test realpath(GAP.Packages.locate_package("hap"))==realpath(selected)
    # Independent fresh Magma/VMH GL complexes, restricted to SL₂×Z, give
    # these homology groups via Künneth. Expected answers live only in tests.
    cases=((D=5,construct=()->hap_voronoi_complex(5;group=:SL,resolution_length=4),
            homology=(Int[],[0,0,0],[2,60,120])),
           (D=12,construct=()->GAP.evalstr("VoronoiComplexSL(2,3 : backend := \"julia\")"),
            homology=([2,6],[0,0,0],[2,2,6,12,24])))
    withenv("PERFECT_FORMS_BACKEND_CONFIG"=>"/no/such/backend/config.json") do
        for c in cases
            K=c.construct()
            @test GAP.Globals.IsHapNonFreeResolution(K)
            metadata=GAP.evalstr("K -> K!.arithmeticVoronoi")(K)
            @test Int(metadata.discriminant)==c.D
            @test Int(metadata.rank)==2
            @test metadata.groupKind==GAP.GapObj("SL")
            @test metadata.levelModel==GAP.GapObj("sections")
            P=arithmetic_problem(c.D)
            generators=GAP.Globals.GeneratorsOfGroup(GAP.evalstr("K -> K!.group")(K))
            @test all(generators) do g
                rows=GAP.gap_to_julia(g;recursive=true)
                A=matrix(QQ,4,4,[QQ(x) for row in rows for x in row])
                det(decode_matrix(P,A))==1
            end
            R=GAP.Globals.FreeGResolution(K,4)
            C=GAP.Globals.TensorWithIntegers(R)
            for (n,expected) in enumerate(c.homology)
                @test sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(C,n))))==expected
            end
        end
    end
    # A requested finite-degree resolution can use a skeleton, but it must
    # never silently advertise missing higher cells as zero chain modules.
    K=hap_voronoi_complex(5;group=:SL,resolution_length=2)
    metadata=GAP.evalstr("K -> K!.arithmeticVoronoi")(K)
    @test metadata.levelModel==GAP.GapObj("sections")
    @test Int(metadata.computedCellDimension)==2
    @test Int(metadata.certifiedCellDimension)==3
    R=GAP.Globals.FreeGResolution(K,2)
    @test isempty(GAP.gap_to_julia(GAP.Globals.Homology(GAP.Globals.TensorWithIntegers(R),1)))
    @test_throws GAP.GAPError GAP.evalstr("K -> K!.dimension(3)")(K)
    @test_throws ArgumentError hap_voronoi_complex(5;group=:SL,rank=1)
    @test_throws ArgumentError hap_voronoi_complex(5;group=:PSL)
    @test_throws ArgumentError hap_voronoi_complex(5;group=:SL,resolution_length=0)
end

@testset "Regular real-quadratic SL incidence through the GAP dispatcher" begin
    # Discriminant 12: an unrestricted face/edge pivot can put four distinct
    # faces around an edge of a single 3-cell. Reconstruct lifts using GAP
    # matrix groups independently of Julia's canonical transporter keys.
    enable_hap_julia_backend(;resolution_length=3)
    K = GAP.evalstr("""VoronoiComplexSL(2,3 : backend := "julia", reduce := "regular")""")
    metadata = GAP.evalstr("K -> K!.arithmeticVoronoi")
    @test metadata(K).reduction == GAP.GapObj("regular")
    two_faces_per_edge = GAP.evalstr("""
        function(K)
            local lift, t, faces, face, edges;
            lift := function(dim, a, g)
                local i, m;
                i := AbsInt(a[1]);
                m := g*K!.elts[a[2]];
                return [i, Minimum(List(Elements(K!.stabilizer(dim,i)), h -> m*h))];
            end;
            for t in [1..K!.dimension(3)] do
                faces := Set(K!.boundary(3,t), a -> lift(2,a,One(K!.group)));
                edges := [];
                for face in faces do
                    Append(edges,Set(K!.boundary(2,face[1]), a -> lift(1,a,face[2])));
                od;
                if IsEmpty(edges) or not ForAll(Collected(edges), p -> p[2]=2) then
                    return false;
                fi;
            od;
            return K!.dimension(3)>0;
        end
    """)
    @test two_faces_per_edge(K)
    R = GAP.Globals.FreeGResolution(K,3)
    C = GAP.Globals.TensorWithIntegers(R)
    @test sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(C,1)))) == [2,6]
    @test Int.(GAP.gap_to_julia(GAP.Globals.Homology(C,2))) == [0,0,0]

    # The general dispatcher and the quadratic convenience constructor share
    # the callback. Pass real GAP booleans and a GAP string through a record call.
    enable_hap_julia_backend(;resolution_length=3, reduce=:regular)
    request = GAP.evalstr("""
        rec(group := "SL", rank := 2,
            field := rec(kind := "quadratic", discriminant := 8),
            order := rec(kind := "maximal"), lattice := rec(kind := "free"))
    """)
    construct = GAP.evalstr("""{r,mode} -> VoronoiComplex(r : backend := "julia", reduce := mode)""")
    orbit_count = GAP.evalstr("K -> Sum([0..3], K!.dimension)")
    original = construct(request,false)
    @test metadata(original).reduction == GAP.GapObj("none")
    for (option,expected) in ((true,"algebraic"), (GAP.GapObj("regular"),"regular"))
        reduced = construct(request,option)
        @test metadata(reduced).reduction == GAP.GapObj(expected)
        @test orbit_count(reduced) < orbit_count(original)
    end

    enable_hap_julia_backend(;resolution_length=2, reduce=:regular)
    partial = GAP.evalstr("""VoronoiComplexSL(2,5 : backend := "julia")""")
    @test metadata(partial).reduction == GAP.GapObj("regular")
    @test Int(metadata(partial).computedCellDimension) == 2
    @test Int(metadata(partial).omittedTopGenerators) == 0
    geometric = GAP.evalstr("""VoronoiComplexSL(2,5 : backend := "julia", reduce := false)""")
    @test metadata(geometric).reduction == GAP.GapObj("none")
    @test_throws ArgumentError enable_hap_julia_backend(;reduce=:unknown)
    @test_throws ArgumentError ArithmeticVoronoi.RealQuadraticSL.compute_real_sl_complex(5;reduce=1)
    enable_hap_julia_backend(;resolution_length=4)
end
