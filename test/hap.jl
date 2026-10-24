using ArithmeticVoronoi, Oscar, Test
import Oscar: GAP

# CI also runs this suite against a pinned upstream checkout. Check its entry
# points before enabling the adapter, which otherwise supports older releases.
hap_test_root=get(ENV,"ARITHMETIC_VORONOI_HAP_ROOT","")
@testset "HAP package selection" begin
    if isempty(hap_test_root)
        @test GAP.Packages.load("hap")
    else
        @test GAP.Packages.load(realpath(hap_test_root))
        @test realpath(GAP.Packages.locate_package("hap"))==realpath(hap_test_root)
        for name in ("RegisterVoronoiBackend","VoronoiComplex","VoronoiComplexGL","VoronoiComplexSL")
            @test GAP.Globals.IsBoundGlobal(GAP.GapObj(name))
        end
    end
end
println("HAP ",GAP.gap_to_julia(GAP.Globals.GAPInfo.PackagesLoaded.hap[2]),
        " from ",GAP.Packages.locate_package("hap"))
flush(stdout)

const HAV=ArithmeticVoronoi
hap_homology(R,n)=sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(
    GAP.Globals.TensorWithIntegers(R),n))))

@testset "native equivariant Voronoi complexes" begin
    # These cell counts and higher homology come from independent fresh VMH
    # computations; the tests themselves enumerate natively and never run Magma.
    expected=Dict(-7=>([1,2,1],[2,2,2],[2,12]),
                   5=>([2,4,5,4,1],[0,2],[2,2,2,2]))
    @test GAP.Packages.load("hap")
    native=Dict()
    hap=Dict()
    withenv("PERFECT_FORMS_BACKEND_CONFIG"=>"/no/such/backend/config.json") do
        for D in (-7,5)
            E=enumerate_perfect_forms(arithmetic_problem(D);backend=:native)
            C=voronoi_complex(E)
            native[D]=C
            @test length.(C.cells)==expected[D][1]
            @test validate_voronoi_complex(C)
            K=hap_complex(C;resolution_length=4)
            hap[D]=K
            R=GAP.Globals.FreeGResolution(K,3)
            @test hap_homology(R,1)==expected[D][2]
            @test hap_homology(R,2)==expected[D][3]
        end
    end
    # Constructing another complex must not overwrite captured data in the first.
    @test hap_homology(GAP.Globals.FreeGResolution(hap[-7],4),3)==[2,2,2,2,12]
    # A changed orientation must be detected before the object can reach HAP.
    C=native[-7]
    broken=deepcopy(C.boundaries)
    j,g=broken[3][1][1]
    broken[3][1][1]=(-j,g)
    bad=VoronoiComplex(C.enumeration,C.group,C.cells,C.stabilizers,broken,C.generators)
    @test_throws ArithmeticValidationError validate_voronoi_complex(bad)
    @test_throws ArithmeticValidationError hap_complex(bad)
    omitted=deepcopy(C.boundaries)
    pop!(omitted[3][1])
    incomplete=VoronoiComplex(C.enumeration,C.group,C.cells,C.stabilizers,omitted,C.generators)
    @test_throws ArithmeticValidationError validate_voronoi_complex(incomplete)
    # The actual determinant-one subgroup, including central -I.
    S=voronoi_complex(C.enumeration;group=:SL)
    @test all(g->det(g)==1,S.generators)
    @test all(layer->all(stab->all(g->det(g)==1,stab),layer),S.stabilizers)
    R=GAP.Globals.FreeGResolution(hap_complex(S;resolution_length=3),3)
    reference=GAP.Globals.ResolutionSL2QuadraticIntegers(-7,3)
    @test hap_homology(R,1)==[0,4]
    @test hap_homology(R,2)==hap_homology(reference,2)
    @test_throws ArgumentError voronoi_complex(native[5].enumeration;group=:SL)
    @test_throws ArgumentError hap_complex(C;resolution_length=0)
    @test_throws ArgumentError hap_voronoi_complex(-12)
end

@testset "imaginary units and nontrivial class groups" begin
    for (D,gl,sl) in ((-3,[6],[3]),(-4,[2,4],[2,2]),
                      (-15,[2,2,2,2],[0,0,12]),(-23,[2,2,2,2,2],[0,0,0,12]))
        E=enumerate_perfect_forms(arithmetic_problem(D);backend=:native)
        for (group,expected) in ((:GL,gl),(:SL,sl))
            V=voronoi_complex(E;group)
            @test validate_voronoi_complex(V)
            R=GAP.Globals.FreeGResolution(hap_complex(V;resolution_length=2),2)
            @test hap_homology(R,1)==expected
        end
    end
end

@testset "optional HAP Julia provider" begin
    entry_points=Dict()
    if GAP.Globals.IsBoundGlobal(GAP.GapObj("RegisterVoronoiBackend"))
        for name in (:RegisterVoronoiBackend,:VoronoiComplex,:VoronoiComplexGL,:VoronoiComplexSL)
            entry_points[name]=getproperty(GAP.Globals,name)
        end
    end
    enable_hap_julia_backend()
    for (name,entry_point) in entry_points
        @test GAP.Globals.IsIdenticalObj(entry_point,getproperty(GAP.Globals,name))
    end
    construct=GAP.evalstr("r -> VoronoiComplex(r : backend := \"julia\")")
    request=voronoi_request(-7)
    Kgeneral=construct(GAP.GapObj(request;recursive=true))
    @test hap_homology(GAP.Globals.FreeGResolution(Kgeneral,3),2)==[2,12]
    Kdirect=hap_voronoi_complex(voronoi_request(-7;group=:SL);resolution_length=3)
    @test hap_homology(GAP.Globals.FreeGResolution(Kdirect,3),2)==[12]
    @test_throws GAP.GAPError GAP.evalstr("VoronoiComplex(rec())")
    @test_throws GAP.GAPError GAP.evalstr("VoronoiComplex(rec() : backend := \"absent\")")
    @test_throws GAP.GAPError GAP.evalstr("VoronoiComplex(2 : backend := \"julia\")")
    @test_throws GAP.GAPError GAP.evalstr("RegisterVoronoiBackend(\"magma\", r -> r)")
    for (key,value) in (("rank",1),("order",Dict("kind"=>"nonmaximal")),
                       ("lattice",Dict("kind"=>"projective")),
                       ("field",Dict("kind"=>"number_field","polynomial"=>[-1,-1,0,1])))
        bad=deepcopy(request);bad[key]=value
        @test_throws GAP.GAPError construct(GAP.GapObj(bad;recursive=true))
    end
    # GAP's interface takes radicands: -1 must request discriminant -4.
    K=GAP.evalstr("VoronoiComplexGL(2,-1 : backend := \"julia\")")
    @test GAP.evalstr("K -> K!.arithmeticVoronoi.discriminant")(K)==-4
    @test hap_homology(GAP.Globals.FreeGResolution(K,2),1)==[2,4]
    enable_hap_julia_backend() # Re-registering/re-reading is safe.
    S=GAP.evalstr("VoronoiComplexSL(2,-7 : backend := \"julia\")")
    @test hap_homology(GAP.Globals.FreeGResolution(S,2),1)==[0,4]
    @test GAP.evalstr("K -> K!.arithmeticVoronoi.groupKind")(S)==GAP.GapObj("SL")
    Kreal=GAP.evalstr("VoronoiComplexGL(2,2 : backend := \"julia\")")
    @test GAP.evalstr("K -> K!.arithmeticVoronoi.discriminant")(Kreal)==8
    @test hap_homology(GAP.Globals.FreeGResolution(Kreal,2),1)==
          hap_homology(GAP.Globals.ResolutionGL2QuadraticIntegers(2,2),1)
    @test_throws GAP.GAPError GAP.evalstr("VoronoiComplexGL(1,-7 : backend := \"julia\")")
    @test_throws GAP.GAPError GAP.evalstr("VoronoiComplexGL(2,-12 : backend := \"julia\")")
    @test_throws GAP.GAPError GAP.evalstr("VoronoiComplexSL(2,5 : backend := \"julia\")")
end

include("higher_rank.jl")
