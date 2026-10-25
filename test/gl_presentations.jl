using ArithmeticVoronoi, Oscar, JSON3, Test
import Oscar: GAP
const GLAV=ArithmeticVoronoi

function test_gl_presentations(; with_hap=false)
    @testset "GL presentations in arbitrary rank" begin
        if with_hap
            @test GAP.Packages.load(get(ENV,"ARITHMETIC_VORONOI_HAP_ROOT","hap"))
        end
        for (D,n,expected) in ((1,2,[2,2]),(1,3,[2]),(-4,3,[4]))
            mktempdir() do directory
                input=D==1 ? classical_problem(;rank=n) : voronoi_request(D;rank=n)
                R=withenv("PERFECT_FORMS_BACKEND_CONFIG"=>"/no/such/backend/config.json") do
                    run_gl_pipeline(input;directory,progress=false)
                end
                P=R.enumeration.problem;S=R.simplified[:GL];I=identity_matrix(P.field,n)
                @test R.run["status"]=="complete"
                @test verify_presentation(R.gl)
                @test R.run["presentation"]["abelianization"].free_rank==0
                @test parse.(Int,R.run["presentation"]["abelianization"].torsion)==expected
                if n==2
                    @test R.run["geometry"]["local_ridges"]==[0]
                    @test R.run["geometry"]["oriented_loops"]==0
                else
                    @test all(r->r.field_rank==n && length(r.rank_witness)==n,Iterators.flatten(R.context.ridges))
                    @test all(g->isnothing(g.multiplication_table),R.gl.vertex_groups)
                    saved=JSON3.read(read(joinpath(directory,"GL_raw.json"),String))
                    @test all(v->haskey(v,:generator_products),saved.vertex_groups)
                    coverage=JSON3.read(read(joinpath(directory,"two_cells/coverage.json"),String),Dict{String,Any})
                    flags=R.gl.ridge_coverage.flags
                    @test coverage["coverage_mode"]=="finite_stabilizer_orbits"
                    @test GLAV.validate_ridge_coverage(R.context,flags,coverage)
                    @test length(flags)<2sum(length,R.context.ridges)
                    @test_throws ArithmeticValidationError GLAV.validate_ridge_coverage(R.context,flags[1:end-1],coverage)
                    bad=deepcopy(coverage);empty!(bad["local_ridge_orbits"][1][1]["members"])
                    @test_throws ArithmeticValidationError GLAV.validate_ridge_coverage(R.context,flags,bad)
                    bad=deepcopy(coverage)
                    # Choose a member whose recorded transporter is nonidentity.
                    members=Iterators.flatten(o["members"] for o in bad["local_ridge_orbits"][1])
                    member=first(m for m in members if GLAV.field_matrix_wire(P,m["to_representative"])!=I)
                    member["stabilizer_word"]=Int[]
                    @test_throws ArithmeticValidationError GLAV.validate_ridge_coverage(R.context,flags,bad)
                    for (i,id) in unique([(i,id) for (i,id,_) in flags])
                        ridge=R.context.ridges[i][id];rep=R.enumeration.representatives[i]
                        oracle=matrix(QQ,rays(positive_hull(rep.evaluation_rays[ridge.indices,:])))
                        @test Set(GLAV.primitive_ray(v) for v in GLAV.matrix_rows(oracle))==
                            Set(GLAV.primitive_ray(v) for v in GLAV.matrix_rows(ridge.rays))
                    end
                end
                elementary=deepcopy(I);elementary[1,n]=P.degree==1 ? P.field(2) : P.generator
                perm=deepcopy(I);perm[1,1]=perm[2,2]=0;perm[1,2]=perm[2,1]=1
                unit=deepcopy(I);unit[n,n]=P.degree==1 ? -one(P.field) : P.generator
                for g in (elementary,perm,unit,-I,elementary*perm*unit*elementary)
                    for simplified in (false,true)
                        c=matrix_to_word(R,g;simplified)
                        @test verify_word_certificate(R,g,c)
                        bad=deepcopy(c);bad["raw_GL_path"]["objective_denominator"]="0"
                        @test_throws ArithmeticValidationError verify_word_certificate(R,g,bad)
                    end
                end
                @test_throws WordReductionError matrix_to_word(R,I;target=:SL)
                @test_throws WordReductionError matrix_to_word(R,I/2)
                if with_hap && n==3
                    # A second construction of the presentation, via HAP's free
                    # resolution of the equivariant complex, uses different generators.
                    V=voronoi_complex(R.enumeration)
                    F=GAP.Globals.FreeGResolution(hap_complex(V;resolution_length=2),2)
                    T=GAP.Globals.PresentationOfResolution(F)
                    data=GAP.gap_to_julia(GAP.evalstr(
                        "{R,P} -> [List(P.gens, i -> R!.elts[i]), List(P.relators, LetterRepAssocWord)]")(F,T);recursive=true)
                    m=n*P.degree
                    matrices=[decode_matrix(P,matrix(QQ,m,m,[A[i][j] for i in 1:m for j in 1:m])) for A in data[1]]
                    relators=[Int.(w) for w in data[2]]
                    @test all(w->GLAV.evaluate_word(matrices,w)==I,relators)
                    for g in matrices
                        @test verify_word_certificate(R,g,matrix_to_word(R,g))
                    end
                    @test parse.(Int,abelianization_invariants(length(matrices),relators).torsion)==expected
                    @test sort(Int.(GAP.gap_to_julia(GAP.Globals.Homology(
                        GAP.Globals.TensorWithIntegers(F),1))))==expected
                end
            end
        end
        @test_throws ArgumentError run_gl_pipeline(voronoi_request(-4;rank=3,group=:SL))
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    test_gl_presentations()
end
