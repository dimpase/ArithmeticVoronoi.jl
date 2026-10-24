using ArithmeticVoronoi, Oscar, JSON3, Test
const TCAV=ArithmeticVoronoi
const TCROOT=normpath(joinpath(@__DIR__,".."))
const TWO_CELL_OBSERVATIONS=Dict{String,Any}()

@testset "exact two-cells, group-action presentations and SL extraction" begin
    for D in (5,8)
        P=arithmetic_problem(D);I=identity_matrix(P.field,2)
        E=enumerate_perfect_forms(P);ctx=TwoCellContext(E)
        directory=joinpath(TCROOT,"results/two_cells","D$D")
        if !TCAV.compatible_two_cell_cache(ctx,directory)
            # Frozen certificates refer to exact generator and facet labels.
            # Another valid enumeration needs its own independently checked loops.
            println("D$D: computing two-cell fixtures for source ",ctx.source_hash);flush(stdout)
            directory=mktempdir(joinpath(TCROOT,"results");prefix="two_cell_fixtures_D$(D)_",cleanup=false)
            enumerate_two_cells(ctx;directory)
        end
        coverage=JSON3.read(read(joinpath(directory,"coverage.json"),String),Dict{String,Any})
        @test coverage["status"]=="complete"
        @test length.(ctx.ridges)==(D==5 ? [75,15] : [15,76]) # measured regression fixtures
        @test all(r.field_rank==2 && length(r.facets)==2 && rank(r.rays)==4 for list in ctx.ridges for r in list)
        certificates=TwoCellCertificate[]
        for entry in coverage["entries"]
            @test entry["status"]=="complete"
            pair=[TCAV.read_two_cell(ctx,joinpath(directory,file)) for file in entry["certificates"]]
            @test length(pair)==2
            @test length(TCAV.reverse_cycle_check(ctx,pair...))==length(first(pair).steps)
            append!(certificates,pair)
        end
        @test length(certificates)==2sum(length,ctx.ridges)
        @test Set((c.source,c.ridge) for c in certificates)==Set((i,r.id) for (i,list) in enumerate(ctx.ridges) for r in list)
        @test all(validate_two_cell(ctx,c) for c in certificates)
        @test any(c.residual==-I for c in certificates)
        first_cell=first(certificates)
        @test length(first_cell.steps)==(D==5 ? 4 : 5)
        # A noncommuting crossing distinguishes b*s*a from a*s*b as actual forms.
        @test any(TCAV.act_form(E.representatives[s.next.representative].form,s.b*s.s*s.state.frame)!=
            TCAV.act_form(E.representatives[s.next.representative].form,s.state.frame*s.s*s.b) for c in certificates for s in c.steps)
        bad_steps=copy(first_cell.steps)
        bad_steps[1]=merge(bad_steps[1],(;b=-bad_steps[1].b))
        bad=TwoCellCertificate(first_cell.source,first_cell.ridge,first_cell.initial_outgoing,bad_steps,
            first_cell.residual,first_cell.residual_word,first_cell.attaching_word,first_cell.source_hash)
        @test_throws TCAV.ArithmeticValidationError validate_two_cell(ctx,bad)
        central=only(findall(g->g==-I,E.representatives[first_cell.source].stabilizer))
        bad=TwoCellCertificate(first_cell.source,first_cell.ridge,first_cell.initial_outgoing,first_cell.steps,
            first_cell.residual,[first_cell.residual_word;central],first_cell.attaching_word,first_cell.source_hash)
        @test_throws TCAV.ArithmeticValidationError validate_two_cell(ctx,bad)
        @test -I!=I
        @test all(TCAV.signed_key(-v)==TCAV.signed_key(v) for v in E.representatives[1].minimal_vectors)

        scratch=mktempdir(joinpath(TCROOT,"results");prefix="two_cell_tests_",cleanup=false)
        data=JSON3.read(JSON3.write(two_cell_data(ctx,first_cell)),Dict{String,Any})
        for (name,mutate) in [("wrong_source",d->(d["source_hash"]="0"^64)),
                               ("incomplete",d->(d["status"]="partial")),
                               ("wrong_group",d->(d["group"]="PGL(2,O)")),
                               ("bad_transporter",d->(d["steps"][1]["b"]=TCAV.field_matrix_data(P,-first_cell.steps[1].b)))]
            modified=deepcopy(data);mutate(modified);path=joinpath(scratch,name*".json");TCAV.atomic_json(path,modified)
            @test_throws TCAV.ArithmeticValidationError TCAV.read_two_cell(ctx,path)
        end
        partial_path=joinpath(scratch,"limited_walk.json")
        @test_throws TCAV.ArithmeticValidationError walk_ridge(ctx,1,1;max_steps=1,checkpoint=partial_path)
        @test JSON3.read(read(partial_path,String))["status"]!="complete"
        if D==5
            batch=enumerate_two_cells(ctx;directory=joinpath(scratch,"resume"),stop_after=1)
            @test batch.ledger["status"]=="partial"
            resumed=enumerate_two_cells(ctx;directory=joinpath(scratch,"resume"),stop_after=2,resume=true)
            @test resumed.ledger["status"]=="partial" && resumed.ledger["cache_hits"]==2
            @test count(e->e["status"]=="complete",resumed.ledger["entries"])==2
            @test any(e->e["status"]=="pending",resumed.ledger["entries"])
        end

        invalid_coverage=deepcopy(coverage);invalid_coverage["entries"][1]["status"]="pending"
        @test_throws TCAV.ArithmeticValidationError build_gl_presentation(ctx,certificates;coverage=invalid_coverage)
        invalid_coverage=deepcopy(coverage);invalid_coverage["source_hash"]="0"^64
        @test_throws TCAV.ArithmeticValidationError build_gl_presentation(ctx,certificates;coverage=invalid_coverage)
        @test_throws TCAV.ArithmeticValidationError build_gl_presentation(ctx,certificates[1:end-2];coverage)
        C=build_gl_presentation(ctx,certificates;coverage)
        @test verify_presentation(C)
        @test all(TCAV.evaluate_word(C.generator_matrices,w)==I for w in C.relators)
        @test all(C.coherent.labels[k]==I for k in C.coherent.tree)
        @test all(-I in g.generators for g in C.vertex_groups)
        @test length(C.graph.edges)==4
        for d in C.graph.darts,u in d.oriented_stabilizer
            @test TCAV.act_form(E.representatives[d.source].form,u)==E.representatives[d.source].form
            v=d.transporter*u*inv(d.transporter)
            @test TCAV.act_form(E.representatives[d.target].form,v)==E.representatives[d.target].form
        end
        if D==8
            inverted=only(e for e in C.graph.edges if !isempty(e.reversals))
            @test length(inverted.oriented_stabilizer)==2 && length(inverted.setwise_stabilizer)==4
            @test all(g*g in inverted.oriented_stabilizer for g in inverted.reversals)
            @test all(g*g==I for g in inverted.reversals) # measured D8 reversals, distinct from the group order
            quarter_turn=matrix(P.field,[0 -1;1 0])
            @test quarter_turn^2==-I && quarter_turn^4==I
            cyclic_four=TCAV.finite_group_record(P,[I,quarter_turn,-I,-quarter_turn])
            @test length(cyclic_four.elements)==4 && -I in cyclic_four.generators
            q=cyclic_four.keys[TCAV.matrix_key(encode_matrix(P,quarter_turn))]
            minus=cyclic_four.keys[TCAV.matrix_key(encode_matrix(P,-I))]
            @test cyclic_four.multiplication_table[q,q]==minus
        else
            @test all(isempty(e.reversals) for e in C.graph.edges)
        end
        # Both search and corrected relations must detect matrix labels, not unsigned permutations.
        constructed=matrix_to_word(C,-I)
        @test TCAV.evaluate_word(C.generator_matrices,constructed.word)==-I
        S=extract_sl_presentation(C)
        @test verify_sl_presentation(S)
        @test S.scalar_reduction.method=="exact separating-facet cone descent"
        @test all(TCAV.rational_wire(s.after)<TCAV.rational_wire(s.before) for s in S.scalar_reduction.steps)
        @test TCAV.evaluate_word(C.generator_matrices,S.scalar_word)==S.units.epsilon*I
        @test length(S.transversal_words)==4
        @test all(det(g)==1 for g in S.generator_matrices)
        @test all(TCAV.evaluate_word(S.generator_matrices,w)==I for w in S.relators)
        @test TCAV.evaluate_word(S.generator_matrices,S.minus_identity_word)==-I
        @test_throws TCAV.ArithmeticValidationError matrix_to_word(C,S.units.epsilon*I;max_cones=1)
        for (label,matrices,relations) in [("GL",C.generator_matrices,C.relators),("SL",S.generator_matrices,S.relators)]
            simple=simplify_matrix_presentation(matrices,relations;label)
            @test verify_simplification(matrices,relations,simple)
            @test simple.isomorphism_verified
            @test length(simple.generator_matrices)<length(matrices)
            bad_maps=deepcopy(simple.forward_generator_words);bad_maps[1]=Int[]
            @test_throws TCAV.ArithmeticValidationError verify_simplification(matrices,relations,merge(simple,(;forward_generator_words=bad_maps)))
        end
        TWO_CELL_OBSERVATIONS[string(D)]=Dict("status"=>"passed","ridges"=>length.(ctx.ridges),"loops"=>length(certificates),
            "GL_generators"=>length(C.generator_matrices),"GL_relators"=>length(C.relators),
            "SL_generators"=>length(S.generator_matrices),"SL_relators"=>length(S.relators),"diagnostics"=>scratch)
    end
end
TCAV.atomic_json(joinpath(TCROOT,"results/two_cell_test_summary.json"),Dict("status"=>"passed","examples"=>TWO_CELL_OBSERVATIONS))
