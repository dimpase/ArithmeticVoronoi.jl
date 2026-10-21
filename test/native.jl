using ArithmeticVoronoi, Oscar, JSON3, Test
const NAV=ArithmeticVoronoi

@testset "exact native deformation" begin
    A=identity_matrix(QQ,2)
    R=matrix(QQ,[0 1;1 0])
    next=NAV.native_flip(A,R)
    @test next.gram==matrix(QQ,[1 1//2;1//2 1])
    @test next.data.minimum==1
    @test length(next.data.vectors)==3
    @test_throws ArithmeticValidationError NAV.native_flip(A,zero_matrix(QQ,2,2))
    @test_throws ArithmeticValidationError NAV.native_flip(A,identity_matrix(QQ,2))
    # Preserving a minimum face is insufficient if the direction is PSD:
    # that ray never acquires a new short vector. Its negative reaches diag(1,1).
    B=matrix(QQ,[1 0;0 2])
    S=matrix(QQ,[0 0;0 1])
    @test_throws ArithmeticValidationError NAV.native_flip(B,S)
    @test NAV.native_flip(B,-S).gram==A
    @test arithmetic_problem(33).discriminant==33
    for D in (-12,0,1,4,9,16,20,45)
        @test_throws ArgumentError arithmetic_problem(D)
    end
    # The trace identity is not perfect; construct the full six-dimensional seed.
    seed=NAV.native_initial_form(arithmetic_problem(8))
    @test seed.data.minimum==1
    @test rank(evaluation_rays(arithmetic_problem(8),seed.data.vectors))==6
end

@testset "native closed neighbour graphs without an executable" begin
    withenv("PERFECT_FORMS_BACKEND_CONFIG"=>"/no/such/backend/config.json") do
        for (D,count) in ((1,1),(5,2),(8,2),(12,3),(13,9))
            P=D==1 ? classical_problem() : arithmetic_problem(D)
            E=enumerate_perfect_forms(P;backend=:native)
            @test E.completion_status==:complete && isempty(E.pending_facets)
            @test length(E.representatives)==count
            @test E.provenance["backend"]=="native"
            @test !E.provenance["from_cache"]
            @test all(c["unresolved_facets"]==0 for c in E.provenance["closure"])
            for i in eachindex(E.representatives),j in 1:i-1
                @test isnothing(equivalence(P,E.representatives[i].form,E.representatives[j].form))
            end
            for r in E.representatives
                @test shortest_data(r.trace_gram).minimum==1
                @test rank(r.evaluation_rays)==length(P.form_basis)
                @test -identity_matrix(P.field,P.rank) in r.stabilizer
                @test finite_stabilizer(P,r.form).order==r.stabilizer_order
                for f in r.facets
                    g=identity_matrix(P.field,P.rank)
                    for k in f.stabilizer_word
                        g*=k>0 ? r.stabilizer[k] : inv(r.stabilizer[-k])
                    end
                    @test g==f.to_orbit_representative
                    @test transpose(g)*r.form*g==r.form
                    index=Dict(Tuple(v)=>k for (k,v) in enumerate(r.minimal_vectors))
                    perm=[index[NAV.signed_key(NAV.apply_integral(encode_matrix(P,g),v))] for v in r.minimal_vectors]
                    @test perm==f.induced_minimal_vector_permutation
                end
            end
            for edge in E.neighbours
                R=encode_matrix(P,edge.transporter)
                @test check_group_matrix(P,R)
                @test transpose(R)*E.representatives[edge.target].trace_gram*R==edge.neighbour_gram
                source=E.representatives[edge.source]
                @test findall(v->NAV.quadratic_value(edge.neighbour_gram,v)==1,source.minimal_vectors)==edge.facet
            end
        end
    end
end

@testset "native enumeration through group presentations and words" begin
    mktempdir() do directory
        R=withenv("PERFECT_FORMS_BACKEND_CONFIG"=>"/no/such/backend/config.json") do
            run_field_pipeline(12;backend=:native,directory,reuse_baseline=false)
        end
        @test R.run["backend"]=="native"
        @test length(R.enumeration.representatives)==3
        @test verify_presentation(R.gl)
        @test verify_sl_presentation(R.sl)
        @test R.run["presentations"]["SL"]["abelianization"].free_rank==0
        @test R.run["presentations"]["SL"]["abelianization"].torsion==["2","6"]
        P=R.enumeration.problem
        g=matrix(P.field,[1 P.generator;0 1])
        for target in (:GL,:SL,:PSL)
            certificate=matrix_to_word(R,g;target)
            @test verify_word_certificate(R,g,certificate)
        end
        # A requested stop must persist current coverage even when periodic
        # ledger writes are throttled. Resuming replays the saved loop files.
        partial_directory=joinpath(directory,"partial_ridges")
        partial=enumerate_two_cells(R.context;directory=partial_directory,stop_after=1)
        @test partial.ledger["status"]=="partial"
        resumed=enumerate_two_cells(R.context;directory=partial_directory,stop_after=2)
        saved=JSON3.read(read(joinpath(partial_directory,"coverage.json"),String))
        @test resumed.ledger["cache_hits"]==saved.cache_hits==2
        @test saved.status=="partial" && saved.completed_oriented_loops==4
        @test count(e->e.status=="complete",saved.entries)==2
        @test any(e->e.status=="pending",saved.entries)
    end
end
