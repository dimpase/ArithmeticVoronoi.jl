using ArithmeticVoronoi, Oscar, JSON, Test, SHA
include(joinpath(@__DIR__,"..","experimental","RealQuadraticSL","RealQuadraticSL.jl"))
using .RealQuadraticSL
using .RealQuadraticSL: validate_level_complex, replay_collapse, MorseLevelComplex
const AV=ArithmeticVoronoi
import Oscar:GAP
function test_homology(K,degree)
    R=GAP.Globals.FreeGResolution(K,degree+1)
    C=GAP.Globals.TensorWithIntegers(R)
    [sort(BigInt.(GAP.gap_to_julia(GAP.Globals.Homology(C,n)))) for n in 1:degree]
end
function elementary_divisors(invariants)
    result=BigInt[]
    for d in invariants
        if iszero(d)
            push!(result,0)
        else
            for (p,e) in factor(ZZ(d));push!(result,BigInt(p)^e);end
        end
    end
    sort(result)
end
function square_determinant_restriction(K,P,U;representation=:native)
    if !GAP.Globals.IsBoundGlobal(GAP.GapObj("AVMorseTestMembership"))
        GAP.evalstr("""
        DeclareAttribute("AVMorseTestMembership",IsGroup);
        InstallMethod(\\in,[IsMultiplicativeElementWithInverse,IsGroup and HasAVMorseTestMembership],
            function(g,H) return AVMorseTestMembership(H)(g); end);
        """)
    end
    function decode(g)
        rows=GAP.gap_to_julia(g;recursive=true)
        B=matrix(QQ,4,4,[QQ(x) for row in rows for x in row])
        representation==:native && return AV.decode_matrix(P,B)
        # VMH uses the row regular representation on (1,sqrt(d)), not O's basis.
        root=P.root_sum==1 ? 2P.generator-1 : P.generator
        rad=div(P.discriminant,P.root_sum==1 ? 1 : 4)
        @assert all(B[2i,2j]==B[2i-1,2j-1] && B[2i,2j-1]==rad*B[2i-1,2j] for i in 1:2 for j in 1:2)
        matrix(P.field,2,2,[P.field(B[2i-1,2j-1])+B[2i-1,2j]*root for i in 1:2 for j in 1:2])
    end
    function encode(g)
        representation==:native && return AV.encode_matrix(P,g)
        rad=div(P.discriminant,P.root_sum==1 ? 1 : 4)
        B=zero_matrix(QQ,4,4)
        for i in 1:2,j in 1:2
            b=coeff(g[i,j],1)/(P.root_sum==1 ? 2 : 1)
            a=coeff(g[i,j],0)+(P.root_sum==1 ? b : 0)
            B[2i-1,2j-1]=a;B[2i-1,2j]=b;B[2i,2j-1]=rad*b;B[2i,2j]=a
        end
        B
    end
    wire(g)=[[BigInt(numerator(x))//BigInt(denominator(x)) for x in row] for row in AV.matrix_rows(encode(g))]
    residue(g)=AV.unit_coordinates(P,U,det(g)).mod_squares
    T=[diagonal_matrix(P.field,[u,one(P.field)]) for u in (one(P.field),-one(P.field),U.epsilon,-U.epsilon)]
    @test residue.(T)==collect(0:3)
    G=GAP.evalstr("K -> K!.group")(K)
    generators=[decode(g) for g in GAP.Globals.GeneratorsOfGroup(G)]
    @test all(g->decode(GAP.GapObj(wire(g);recursive=true))==g,generators)
    schreier=unique([t*g*inv(T[residue(t*g)+1]) for t in T for g in generators])
    @test all(g->residue(g)==0,schreier)
    H=GAP.Globals.Group(GAP.GapObj(wire.(schreier);recursive=true))
    GAP.Globals.SetAVMorseTestMembership(H,GAP.GapObj(g->residue(decode(g))==0))
    result=GAP.Globals.RestrictNonFreeResolution(K,H)
    @test GAP.evalstr("K -> Length(K!.cosetRepresentatives)")(result)==4
    result
end

function validate_comparisons(obj;directory,degree=3,magma_directory=nothing)
    P=obj.M.complex.enumeration.problem;D=P.discriminant
    receipt=Dict{String,Any}("D"=>D,"actual_SL_complex"=>true,
        "cell_orbits"=>length.(obj.L.cells),"reduced_module_orbits"=>length.(obj.reduced.cells),
        "certified_links"=>length(obj.certificate.reports),"source"=>"fresh native enumeration; no precomputed answers")
    @testset "Real SL2 D=$D" begin
        @test validate_level_complex(obj.L)
        @test validate_level_complex(obj.reduced)
        @test length(obj.certificate.reports)==2sum(length,obj.M.complex.cells)
        # Reject an incomplete contraction witness and a corrupted orientation.
        r=first(x for x in obj.certificate.reports if !isempty(x.result.pairs))
        link=obj.data["links"][r.orbit]
        delta=[obj.certificate.bases[v["orbit"]]+v["shift"]-obj.certificate.bases[r.orbit] for v in link["vertices"]]
        facets=[Tuple(v for v in f if r.sgn*delta[v]>0) for f in link["facets"]]
        filter!(!isempty,facets)
        @test replay_collapse(facets,r.result.pairs)
        @test_throws AssertionError replay_collapse(facets,r.result.pairs[2:end])
        b=[[copy(x) for x in layer] for layer in obj.L.boundaries]
        j,g=b[2][1][1];b[2][1][1]=(-j,g)
        bad=MorseLevelComplex(obj.L.experiment,obj.L.certificate,obj.L.flag_orbits,obj.L.cells,obj.L.stabilizers,b)
        @test_throws AssertionError validate_level_complex(bad)
        h=test_homology(obj.K,degree)
        receipt["SL_homology"]=h
        println("D$D direct SL homology: ",h);flush(stdout)
        # A separate presentation algorithm computes H1 from relator exponent sums.
        ctx=TwoCellContext(obj.M.complex.enumeration)
        two=enumerate_two_cells(ctx;directory=joinpath(directory,"D$(D)_presentation_checks"),reduce_orbits=false,resume=false)
        pres=extract_sl_presentation(build_gl_presentation(ctx,two.certificates;coverage=two.ledger))
        ab=AV.abelianization_invariants(length(pres.generator_matrices),pres.relators)
        hp=sort([zeros(BigInt,ab.free_rank);parse.(BigInt,ab.torsion)])
        @test elementary_divisors(h[1])==elementary_divisors(hp)
        receipt["presentation_H1"]=hp
        Kgl=hap_complex(obj.M.complex;resolution_length=degree+1)
        H=square_determinant_restriction(Kgl,P,obj.M.units)
        product=test_homology(H,degree)
        receipt["native_SL_times_Z_homology"]=product
        for n in 1:degree
            previous=n==1 ? BigInt[0] : h[n-1]
            @test elementary_divisors(product[n])==sort([elementary_divisors(h[n]);elementary_divisors(previous)])
        end
        if !isnothing(magma_directory)
            root=joinpath(magma_directory,"magma_D$(D)")
            @test isfile(joinpath(root,"complete.txt"))
            GAP.Globals.Read(GAP.GapObj(joinpath(root,"GAP","ReadWellRoundedComplex.gi")))
            Kmag=GAP.Globals.ResolutionFromWellRoundedComplex(GAP.GapObj(joinpath(root,"complex_gl.g")),degree+1)
            Hmag=square_determinant_restriction(Kmag,P,obj.M.units;representation=:magma)
            reference=test_homology(Hmag,degree)
            @test reference==product
            receipt["Magma_SL_times_Z_homology"]=reference
            receipt["Magma_GL_complex_sha256"]=bytes2hex(sha256(read(joinpath(root,"complex_gl.g"))))
            receipt["Magma_reference_directory"]=root
        end
        receipt["status"]="passed"
        write(joinpath(directory,"D$(D)_validated_receipt.json"),JSON.json(receipt;pretty=4))
    end
    receipt
end

function run_real_sl_tests(directory;magma_directory=nothing,discriminants=(5,12))
    mkpath(directory)
    for D in discriminants
        obj=compute_real_sl_complex(D;progress=true,certificate_directory=joinpath(directory,"D$(D)_construction"))
        validate_comparisons(obj;directory,magma_directory)
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (1,2) || error("provide a durable output directory and optionally the fresh Magma reference directory")
    run_real_sl_tests(abspath(ARGS[1]);magma_directory=length(ARGS)==2 ? abspath(ARGS[2]) : nothing)
end
