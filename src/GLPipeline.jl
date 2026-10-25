export GLGroupResult, run_gl_pipeline

struct GLGroupResult
    enumeration::ArithmeticVoronoiEnumeration
    context::TwoCellContext
    gl::PresentationCertificate
    simplified::Dict{Symbol,Any}
    metadata::Any
    directory::String
    run::Dict{String,Any}
end

"""
    run_gl_pipeline(request_or_problem; directory=nothing, simplify=true, progress=true, resume=false)

Compute a certified GL presentation in any supported rank, with exact ridge
cycles, generator matrices and optional Tietze maps. Quadratic request records
and rational `classical_problem(; rank=n)` models use native enumeration.
The separate GL/SL/PSL `run_field_pipeline` retains its rank-two scope.
Saved attaching loops are read only with `resume=true`, then exactly replayed.
"""
function run_gl_pipeline(request::VoronoiRequestRecord;kwargs...)
    run_gl_pipeline(arithmetic_problem(enumeration_request(request));kwargs...)
end
run_gl_pipeline(D::Integer;rank=2,kwargs...)=run_gl_pipeline(arithmetic_problem(Int(D);rank);kwargs...)

function run_gl_pipeline(P::ArithmeticVoronoiProblem;directory=nothing,simplify=true,progress=true,resume::Bool=false)
    label="GL$(P.rank)_D$(P.discriminant)"
    isnothing(directory) && (directory=joinpath(pwd(),"results/gl_presentations",label))
    directory=abspath(directory);mkpath(directory)
    started=time_ns();timings=Dict{String,Any}()
    run=Dict{String,Any}("schema"=>"av_gl_run_v1","status"=>"running","problem"=>problem_descriptor(P),
        "backend"=>"native","stage_seconds"=>timings,"julia_version"=>string(VERSION),
        "oscar_version"=>string(pkgversion(Oscar)),"started_utc"=>string(now(UTC)),"resume_requested"=>resume)
    save()=atomic_json(joinpath(directory,"run.json"),run)
    function stage(name,f)
        run["current_stage"]=name;save()
        println(label,": ",name);flush(stdout);t=time_ns()
        try
            value=f();timings[name]=(time_ns()-t)/1e9;save();value
        catch err
            timings[name]=(time_ns()-t)/1e9;run["status"]="incomplete"
            run["error"]=sprint(showerror,err);save();rethrow()
        end
    end
    metadata=problem_descriptor(P);atomic_json(joinpath(directory,"field.json"),metadata)
    E=stage("enumeration",()->enumerate_perfect_forms(P;backend=:native,progress))
    write_enumeration(joinpath(directory,"enumeration.json"),E)
    ctx=stage("ridges",()->TwoCellContext(E))
    cells=stage("attaching_loops",()->enumerate_two_cells(ctx;directory=joinpath(directory,"two_cells"),reduce_orbits=P.rank>2,resume))
    run["replayed_loop_certificates"]=get(cells.ledger,"cache_hits",0)
    C=stage("GL_assembly",()->build_gl_presentation(ctx,cells.certificates;coverage=cells.ledger))
    atomic_json(joinpath(directory,"GL_raw.json"),presentation_data(C))
    small=Dict{Symbol,Any}()
    if simplify
        small[:GL]=stage("GL_simplification",()->simplify_matrix_presentation(C.generator_matrices,C.relators;
            label="GL($(P.rank),O), D=$(P.discriminant)"))
        atomic_json(joinpath(directory,"GL_simplified.json"),simplification_data(P,small[:GL]))
    end
    object=get(small,:GL,C)
    invariants=stage("abelianization",()->abelianization_invariants(length(object.generator_matrices),object.relators))
    run["presentation"]=Dict("raw_generators"=>length(C.generator_matrices),"raw_relators"=>length(C.relators),
        "generators"=>length(object.generator_matrices),"relators"=>length(object.relators),"abelianization"=>invariants)
    run["geometry"]=Dict("perfect_form_orbits"=>length(E.representatives),"local_ridges"=>length.(ctx.ridges),
        "oriented_loops"=>length(cells.certificates),"local_ridge_orbits"=>cells.ledger["locally_reduced_ridge_counts"],
        "global_unoriented_edges"=>length(C.graph.edges))
    run["source_hash"]=ctx.source_hash;run["total_in_process_seconds"]=(time_ns()-started)/1e9
    run["status"]="complete";save()
    GLGroupResult(E,ctx,C,small,metadata,directory,run)
end
