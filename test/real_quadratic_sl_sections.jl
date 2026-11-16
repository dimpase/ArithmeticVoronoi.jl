using ArithmeticVoronoi, Oscar, Test
const SectionRQSL=ArithmeticVoronoi.RealQuadraticSL

function section_test_homology(K,degree)
    GAP=Oscar.GAP
    R=GAP.Globals.FreeGResolution(K,degree+1)
    C=GAP.Globals.TensorWithIntegers(R)
    [sort(BigInt.(GAP.gap_to_julia(GAP.Globals.Homology(C,n)))) for n in 1:degree]
end

@testset "Certified sections and compatible subdivision" begin
    for D in (5,12)
        P=arithmetic_problem(D)
        C=voronoi_complex(enumerate_perfect_forms(P;backend=:native))
        M=SectionRQSL.morse_experiment(C;maximal_flags=false)
        data,_,faces=SectionRQSL.morse_link_graphs(M)
        certificate=SectionRQSL.certify_morse_graphs(data)
        plan=SectionRQSL.certify_level_sections(M,certificate;faces)
        @test all(r->r.result.passed && r.core_dimension==0,plan.proofs)
        @test SectionRQSL.section_parent_dimension(plan,4,2)==3
        flag=SectionRQSL.build_level_complex(M,data,certificate;faces)
        expected=section_test_homology(SectionRQSL.level_hap_complex(SectionRQSL.reduce_level_complex(flag);resolution_length=4),3)
        for mode in (:sections,:mixed,:flags)
            refine=mode==:sections ? Tuple{Int,Int,Int}[] :
                [(r.parent_dimension+1,r.parent_orbit,r.exponent) for r in plan.sections
                    if mode==:flags || (r.parent_dimension==2 && r.parent_orbit==1)]
            L=SectionRQSL.assemble_section_level(M,certificate,plan;faces,top=3,parent_max_dimension=4,extra_refinements=refine)
            @test SectionRQSL.validate_level_complex(L)
            if D==12 && mode==:sections
                # These reversing elements occur in dimension three; a
                # two-dimensional partial complex correctly omits them.
                has_reversing=any(SectionRQSL.level_orientation(L,k,i,g)==-1
                    for (k,layer) in enumerate(L.stabilizers) for (i,S) in enumerate(layer) for g in S)
                @test has_reversing
            end
            mode==:flags && @test length.(L.cells)==length.(flag.cells)
            if mode==:mixed
                @test any(c->c.kind==:section,Iterators.flatten(L.cells))
                @test any(c->c.kind==:flag,Iterators.flatten(L.cells))
            end
            reduced=SectionRQSL.reduce_level_complex(L)
            @test SectionRQSL.validate_level_complex(reduced)
            K=SectionRQSL.level_hap_complex(reduced;resolution_length=4)
            metadata=Oscar.GAP.evalstr("K -> K!.arithmeticVoronoi")(K)
            @test metadata.levelModel==Oscar.GAP.GapObj("sections")
            @test Int(metadata.retainedOriginalDimension)==4
            @test section_test_homology(K,3)==expected
            # The checker must reject a corrupted incidence before HAP sees it.
            entry=L.boundaries[2][1][1]
            L.boundaries[2][1][1]=(-entry[1],entry[2])
            @test_throws AssertionError SectionRQSL.validate_level_complex(L)
            L.boundaries[2][1][1]=entry
        end
        partial=SectionRQSL.assemble_section_level(M,certificate,plan;faces,top=2,parent_max_dimension=3)
        @test SectionRQSL.validate_level_complex(partial)
        K=SectionRQSL.level_hap_complex(SectionRQSL.reduce_level_complex(partial);resolution_length=2)
        @test only(section_test_homology(K,1))==first(expected)
        @test_throws ArgumentError SectionRQSL.level_hap_complex(partial;resolution_length=3)
        @test_throws Oscar.GAP.GAPError Oscar.GAP.evalstr("K -> K!.dimension(3)")(K)
        metadata=Oscar.GAP.evalstr("K -> K!.arithmeticVoronoi")(K)
        @test Int(metadata.retainedOriginalDimension)==3
    end
    @test_throws ArgumentError SectionRQSL.compute_real_sl_complex(5;cell_model=:unknown)
    # An uncertified low-codimension core cannot be omitted from a partial model.
    plan=(proofs=[(parent=(4,1),core_dimension=1),(parent=(5,1),core_dimension=2),
                  (parent=(8,1),core_dimension=2)],)
    @test SectionRQSL.section_parent_dimension(plan,9,2)==5
    @test SectionRQSL.section_parent_dimension(plan,9,4)==5
    @test SectionRQSL.section_parent_dimension(plan,9,5)==8
end

@testset "Cancellation with noncentral conjugate stabilizers" begin
    P=arithmetic_problem(5;rank=3)
    I=identity_matrix(P.field,P.rank)
    r=diagonal_matrix(P.field,[-1,-1,1])
    g=deepcopy(I);g[1,3]=1
    lower=Any[I,r]
    upper=Any[g*h*inv(g) for h in lower]
    @test Set(lower)!=Set(upper)
    M=SectionRQSL.MorseExperiment((enumeration=(problem=P,),),nothing,nothing,
        Dict{Any,SectionRQSL.MorseVertex}(),Dict{Any,SectionRQSL.MorseVertex}(),
        IdDict{Any,Any}(),Dict{Any,Int}(),Vector{BigInt}[],nothing)
    cell=SectionRQSL.LevelCell(0,I,0,1)
    L=SectionRQSL.MorseLevelComplex(M,nothing,SectionRQSL.FlagOrbit[],
        [[cell,cell],[cell]],[[lower,upper],[upper]],
        [[Tuple{Int,Any}[],Tuple{Int,Any}[]],[Tuple{Int,Any}[(1,g),(-2,I)]]])
    @test SectionRQSL.validate_level_complex(L)
    R=SectionRQSL.reduce_level_complex(L)
    @test length.(R.cells)==[1,0]
    @test SectionRQSL.validate_level_complex(R)
end

@testset "Cancellation updates in several chain degrees" begin
    P=arithmetic_problem(5;rank=3);I=identity_matrix(P.field,P.rank)
    M=SectionRQSL.MorseExperiment((enumeration=(problem=P,),),nothing,nothing,
        Dict{Any,SectionRQSL.MorseVertex}(),Dict{Any,SectionRQSL.MorseVertex}(),
        IdDict{Any,Any}(),Dict{Any,Int}(),Vector{BigInt}[],nothing)
    cell=SectionRQSL.LevelCell(0,I,0,1)
    for n in 3:6
        simplices=[[findall(i->!iszero(mask & (1<<(i-1))),1:n)
            for mask in reverse(1:(1<<n)-1) if count_ones(mask)==k] for k in 1:n]
        indices=[Dict(Tuple(s)=>i for (i,s) in enumerate(layer)) for layer in simplices]
        cells=[[cell for _ in layer] for layer in simplices]
        stabs=[[Any[I] for _ in layer] for layer in simplices]
        boundaries=[[k==1 ? Tuple{Int,Any}[] : Tuple{Int,Any}[
            ((-1)^(j-1)*indices[k-1][Tuple(s[l] for l in eachindex(s) if l!=j)],I)
            for j in eachindex(s)] for s in layer] for (k,layer) in enumerate(simplices)]
        for full in (true,false)
            last=full ? n : n-1
            L=SectionRQSL.MorseLevelComplex(M,nothing,SectionRQSL.FlagOrbit[],
                cells[1:last],stabs[1:last],boundaries[1:last])
            @test SectionRQSL.validate_level_complex(L)
            R=SectionRQSL.reduce_level_complex(L)
            expected=zeros(Int,last);expected[1]=1;full || (expected[end]=1)
            @test length.(R.cells)==expected
            @test SectionRQSL.validate_level_complex(R)
        end
    end
end
