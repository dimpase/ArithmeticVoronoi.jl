using ArithmeticVoronoi,Oscar,JSON3,Test
const FWAV=ArithmeticVoronoi

function run_fields_words_tests(results)
    observations=Dict{String,Any}()
    @testset "fields, constructive words, projective signs and literature maps" begin
        for D in sort!(collect(keys(results)))
            R=results[D];P=R.enumeration.problem;K=P.field;w=P.generator;I=identity_matrix(K,2);e=R.sl.units.epsilon
            @test discriminant(P.order)==D==P.root_sum^2-4P.root_product
            @test R.metadata.ordinary_class_number==(D==40 ? "2" : "1")
            @test abs(R.sl.units.epsilon_exponent)==1
            @test norm(e)==(D==12 ? 1 : -1)
            @test all(R.run["status"][s]=="complete" for s in ("enumeration","two_cells","GL","SL","PSL"))
            @test R.run["geometry"]["oriented_loops"]==2sum(length,R.context.ridges)
            U=matrix(K,[1 1;0 1]);W=matrix(K,[1 w;0 1]);L=matrix(K,[1 0;1 1]);Lw=matrix(K,[1 0;w 1]);H=matrix(K,[e 0;0 inv(e)])
            inputs=Any[("identity",I),("minus_identity",-I),("stabilizer",R.enumeration.representatives[1].stabilizer[1]),
                ("scalar",e*I),("inverse_scalar",inv(e)*I),("unit_diagonal",matrix(K,[e 0;0 1])),
                ("H",H),("U",U),("W",W),("lower_1",L),("lower_w",Lw)]
            for k in -2:2
                push!(inputs,("product_$k",(U*inv(Lw))^k*H^(k+1)*W^(2-k)))
            end
            scratch=mktempdir(joinpath(R.directory);prefix="word_tests_",cleanup=false)
            leaving=nothing;count=0
            for (name,g) in inputs,target in (:GL,:SL,:PSL)
                target!=:GL && det(g)!=1 && continue
                path=joinpath(scratch,"$(name)_$target.json")
                data=matrix_to_word(R,g;target,cache_path=path)
                @test data["status"]=="complete" && !data["cache_replay"]
                @test verify_word_certificate(R,g,data)
                matrices=R.simplified[target].generator_matrices
                value=FWAV.evaluate_word(matrices,Int.(data["word"]))
                @test target==:PSL ? value in (g,-g) : value==g
                replay=matrix_to_word(R,g;target,cache_path=path)
                @test replay["cache_replay"] && replay["word"]==data["word"]
                @test all(FWAV.rational_wire(s["after"])<FWAV.rational_wire(s["before"]) for s in data["raw_GL_path"]["steps"])
                if !isempty(data["raw_GL_path"]["steps"]) && target==:GL
                    leaving=(g,data)
                end
                FWAV.atomic_json(joinpath(R.directory,"words","$(name)_$target.json"),data)
                count+=1
            end
            @test !isnothing(leaving)
            for target in (:GL,:SL,:PSL)
                data=matrix_to_word(R,U;target,simplified=false)
                @test verify_word_certificate(R,U,data)
            end
            @test matrix_to_word(R,encode_matrix(P,U);target=:SL,certificate=false)==matrix_to_word(R,U;target=:SL,certificate=false)
            for bad in (2I,matrix(K,[1//2 0;0 2]),zero_matrix(K,2,3))
                err=try matrix_to_word(R,bad;target=:GL);nothing catch caught;caught end
                @test err isa WordReductionError && err.kind==:invalid_input
            end
            semilinear=block_diagonal_matrix([matrix(QQ,[1 P.root_sum;0 -1]),matrix(QQ,[1 P.root_sum;0 -1])])
            err=try matrix_to_word(R,semilinear;target=:GL);nothing catch caught;caught end
            @test err isa WordReductionError && err.kind==:invalid_input
            @test_throws WordReductionError matrix_to_word(R,matrix(K,[e 0;0 1]);target=:SL)
            partial=joinpath(scratch,"exhausted.json")
            err=try matrix_to_word(R,e*I;target=:GL,max_steps=0,cache_path=partial);nothing catch caught;caught end
            @test err isa WordReductionError && err.kind==:resource_exhausted
            @test !isfile(partial) && JSON3.read(read(partial*".partial",String))["status"]=="incomplete"
            err=try matrix_to_word(R,e*I;target=:GL,timeout=0.0);nothing catch caught;caught end
            @test err isa WordReductionError && err.kind==:resource_exhausted
            g,data=leaving;bad=deepcopy(data)
            old=FWAV.field_matrix_wire(P,bad["raw_GL_path"]["steps"][1]["actual_frame"])
            bad["raw_GL_path"]["steps"][1]["actual_frame"]=FWAV.field_matrix_data(P,-old)
            badpath=joinpath(scratch,"corrupt.json");FWAV.atomic_json(badpath,bad)
            err=try matrix_to_word(R,g;target=:GL,cache_path=badpath);nothing catch caught;caught end
            @test err isa WordReductionError && err.kind==:invalid_certificate
            @test_throws FWAV.ArithmeticValidationError verify_word_certificate(R,g,bad)
            bad_residual=deepcopy(data)
            minus=FWAV.finite_element_word(R.gl,1,-I)
            append!(bad_residual["raw_GL_path"]["residual_stabilizer_word"],minus)
            @test_throws FWAV.ArithmeticValidationError verify_word_certificate(R,g,bad_residual)
            @test FWAV.evaluate_word(R.sl.generator_matrices,R.sl.minus_identity_word)==-I
            @test any(FWAV.evaluate_word(R.simplified[:PSL].generator_matrices,r)==-I for r in R.simplified[:PSL].relators)
            # A target reverse-facet correction r and a source expansion s are
            # different data. Check both labelled maps in the coherent lifts.
            for d in R.gl.graph.darts
                i,j=d.source,d.target;li=R.gl.coherent.lifts[i];lj=R.gl.coherent.lifts[j]
                target_indices=sort!(abs.(R.context.neighbours[i][d.local_orbit].facet_map))
                f=only(f for f in R.enumeration.representatives[j].facets if f.indices==target_indices)
                r=f.to_orbit_representative
                @test FWAV.act_form(R.enumeration.representatives[j].form,r)==R.enumeration.representatives[j].form
                @test FWAV.evaluate_word(R.gl.generator_matrices,R.gl.transporter_words[d.id])==inv(lj)*d.transporter*li
                for source_facet in R.enumeration.representatives[i].facets
                    source_facet.neighbour_orbit==d.local_orbit || continue
                    s=source_facet.to_orbit_representative
                    sword=FWAV.finite_element_word(R.gl,i,inv(li)*s*li)
                    @test FWAV.evaluate_word(R.gl.generator_matrices,[R.gl.transporter_words[d.id];sword])==inv(lj)*d.transporter*s*li
                end
            end
            if D in (5,8,12,13)
                comparison=compare_kw(R)
                @test comparison.status=="certified_by_independent_faithfulness"
                ref=kw_reference(P)
                @test all(FWAV.evaluate_word(ref.generator_matrices,r)==I for r in ref.SL_relators)
                @test all(FWAV.evaluate_word(ref.generator_matrices,r) in (I,-I) for r in ref.PSL_relators)
                if D==8
                    # GAP's simplified basis can vary by platform. These fixed
                    # matrices give the audited identities below; certify their
                    # words in whichever presentation this run computed.
                    a=matrix(K,[-w -1;1-w -1]);b=matrix(K,[1 w;1-w w-1])
                    @test verify_word_certificate(R,a,matrix_to_word(R,a;target=:SL))
                    @test verify_word_certificate(R,b,matrix_to_word(R,b;target=:SL))
                    z=b^4;Pt=z*a*b;Q0=z*a^2*b^2*a;Rt=inv(a)*b*inv(a)*inv(b)
                    @test z==-I && Pt==matrix(K,[1 1+w;0 1]) && Q0==matrix(K,[1 2+w;0 1])
                    @test Q0*inv(Pt)==U && Pt*inv(U)==W
                    @test inv(Rt*Q0)==H
                    S=(Rt*Q0)*a*b^2*U
                    @test z*S==ref.generator_matrices[4]
                    erroneous=[3,2,-3,-1,-1,-2,-2]
                    @test FWAV.evaluate_word(ref.generator_matrices,erroneous)==matrix(K,[1 3e;0 1])
                    @test FWAV.evaluate_word(ref.generator_matrices,erroneous) ∉ (I,-I)
                end
                if D in (5,8)
                    expected=D==5 ? String[] : ["2","2"]
                    @test R.run["presentations"]["SL"]["abelianization"].torsion==expected
                    @test R.run["presentations"]["GL"]["abelianization"].free_rank==1
                end
            else
                @test_throws FWAV.ArithmeticValidationError kw_reference(P)
                ideal40=ideal(P.order,2)+ideal(P.order,P.order(w))
                @test norm(ideal40)==2 && !is_principal(ideal40)
                v=BigInt[2,0,0,1]
                @test foldl(gcd,v)==1
                L40=integer_lattice(;gram=trace_gram(P,I))
                @test FWAV.quadratic_value(gram_matrix(L40),v)==28
                @test any(FWAV.signed_key(x)==Tuple(v) for (x,n) in short_vectors(L40,QQ(28)))
                @test all(mod(x^2,5) ∉ (2,3) for x in 0:4)
            end
            observations[string(D)]=Dict("status"=>"pending_testset_completion","distinct_inputs"=>length(inputs),"certified_word_examples"=>count,"diagnostics"=>scratch,
                "fresh_and_replay"=>true,"raw_and_simplified_alphabets"=>true)
        end
    end
    for (D,R) in results
        observations[string(D)]["status"]="passed"
        R.run["status"]["constructive_words"]="validated_GL_SL_PSL"
        FWAV.atomic_json(joinpath(R.directory,"word_validation.json"),observations[string(D)])
        FWAV.atomic_json(joinpath(R.directory,"run.json"),R.run)
    end
    FWAV.atomic_json(joinpath(FWAV.PACKAGE_ROOT,"results/fields_words_tests.json"),Dict("status"=>"passed","examples"=>observations))
end
