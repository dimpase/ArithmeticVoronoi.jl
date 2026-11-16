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
