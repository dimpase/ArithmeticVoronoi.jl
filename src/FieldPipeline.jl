export ArithmeticGroupResult, run_field_pipeline

struct ArithmeticGroupResult
    enumeration::ArithmeticVoronoiEnumeration
    context::TwoCellContext
    gl::PresentationCertificate
    sl::SLPresentationCertificate
    psl::Any
    simplified::Dict{Symbol,Any}
    metadata::Any
    directory::String
    run::Dict{String,Any}
end

"""
    run_field_pipeline(D_or_request_or_problem; rank=2, directory=nothing, ...)

Compute GL, SL and PSL presentations and constructive word maps on the free
lattice O^n. Use a quadratic discriminant/request, or `classical_problem(rank=n)`
for Z. PSL quotients by the entire scalar centre of SL. Enumeration and loop
construction run afresh unless explicit checkpoint replay is requested.
"""
function run_field_pipeline(D::Integer;rank=2,kwargs...)
    run_field_pipeline(voronoi_request(D;rank);kwargs...)
end
function run_field_pipeline(request::VoronoiRequestRecord;kwargs...)
    request=enumeration_request(request)
    run_field_pipeline(arithmetic_problem(request);kwargs...)
end
function run_field_pipeline(P::ArithmeticVoronoiProblem;directory=nothing,
                            simplify=true,backend=:native,timeout=backend==:native ? Inf : 600.0,
                            resume::Bool=false,progress=false)
    D=P.discriminant
    request=P.degree==1 ? nothing : voronoi_request(D;rank=P.rank)
    label=P.rank==2 ? "D$D" : "rank$(P.rank)/D$D"
    isnothing(directory) && (directory=joinpath(pwd(),"results/fields",label))
    mkpath(directory);start=time_ns();timings=Dict{String,Float64}()
    if isfile(joinpath(directory,"run.json"))
        history=joinpath(directory,"history");mkpath(history)
        cp(joinpath(directory,"run.json"),joinpath(history,"run_$(time_ns()).json"))
    end
    status=Dict{String,Any}("enumeration"=>"pending","two_cells"=>"pending","GL"=>"pending",
        "SL"=>"pending","PSL"=>"pending","simplification"=>Dict(k=>(simplify ? "pending" : "not_requested") for k in ("GL","SL","PSL")),"constructive_words"=>"unvalidated",
        "literature_comparison"=>"not_requested")
    run=Dict{String,Any}("schema"=>"av_field_run_v1","discriminant"=>D,"rank"=>P.rank,"request"=>request,"status"=>status,
        "stage_seconds"=>timings,"julia_version"=>string(VERSION),"oscar_version"=>string(pkgversion(Oscar)),
        "julia_threads"=>Threads.nthreads(),"started_utc"=>string(now(UTC)),"resume_requested"=>resume)
    save()=atomic_json(joinpath(directory,"run.json"),run)
    stage(label,f)=begin
        println("D$D: ",label);flush(stdout)
        run["current_stage"]=label;run["current_stage_status"]="running";save()
        t=time_ns()
        try
            value=f();timings[label]=(time_ns()-t)/1e9
            run["current_stage_status"]="complete";save();value
        catch err
            timings[label]=(time_ns()-t)/1e9;run["current_stage_status"]="incomplete"
            run["error"]=sprint(showerror,err);save();rethrow()
        end
    end
    save()
    metadata=stage("field_and_unit_verification",()->field_metadata(P));atomic_json(joinpath(directory,"field.json"),metadata)
    E=stage("perfect_enumeration_and_validation",()->enumerate_perfect_forms(P;backend,progress,timeout,force=!resume))
    run["backend"]=string(backend)
    status["enumeration"]="complete";run["enumeration_from_cache"]=E.provenance["from_cache"]
    run["backend_seconds"]=get(E.provenance,"backend_seconds",nothing)
    write_enumeration(joinpath(directory,"enumeration.json"),E);save()
    ctx=stage("ridges",()->TwoCellContext(E))
    certdir=joinpath(directory,"two_cells")
    cells=stage("two_cells_and_coverage",()->enumerate_two_cells(ctx;directory=certdir,reduce_orbits=P.rank>2,resume))
    run["replayed_loop_certificates"]=get(cells.ledger,"cache_hits",0)
    status["two_cells"]="complete";run["certificate_directory"]=relpath(certdir,directory)
    C=stage("GL_assembly",()->build_gl_presentation(ctx,cells.certificates;coverage=cells.ledger))
    atomic_json(joinpath(directory,"GL_raw.json"),presentation_data(C));status["GL"]="complete";save()
    S=stage("SL_extraction",()->extract_sl_presentation(C;timeout))
    atomic_json(joinpath(directory,"SL_raw.json"),sl_presentation_data(S));status["SL"]="complete";save()
    Q=stage("PSL_quotient",()->psl_presentation(S))
    atomic_json(joinpath(directory,"PSL_raw.json"),psl_presentation_data(P,Q));status["PSL"]="complete";save()
    small=Dict{Symbol,Any}()
    for (name,object) in [(:GL,C),(:SL,S),(:PSL,Q)]
        simplify || continue
        path=joinpath(directory,string(name)*"_simplified.json")
        small[name]=stage(string(name)*"_simplification",()->
            simplify_matrix_presentation(object.generator_matrices,object.relators;
                label="$name($(P.rank),O), D$D",projective=name==:PSL))
        atomic_json(path,simplification_data(P,small[name]));status["simplification"][string(name)]="complete";save()
    end
    run["geometry"]=Dict("perfect_form_orbits"=>length(E.representatives),"stabilizer_orders"=>[r.stabilizer_order for r in E.representatives],
        "local_directed_facet_orbits"=>length(E.neighbours),"global_unoriented_edges"=>length(C.graph.edges),
        "local_ridges"=>length.(ctx.ridges),"locally_reduced_ridges"=>cells.ledger["locally_reduced_ridge_counts"],
        "oriented_loops"=>length(cells.certificates),"inverted_edges"=>count(e->!isempty(e.reversals),C.graph.edges),
        "global_two_cell_orbits"=>nothing,"global_two_cell_reduction"=>"not performed; exhaustive redundant local coverage")
    run["presentations"]=stage("abelianizations",()->Dict(string(name)=>Dict("raw_generators"=>length(object.generator_matrices),"raw_relators"=>length(object.relators),
        "simplified_generators"=>haskey(small,name) ? length(small[name].generator_matrices) : nothing,
        "simplified_relators"=>haskey(small,name) ? length(small[name].relators) : nothing,
        "abelianization_presentation"=>haskey(small,name) ? "simplified through the certified Tietze maps" : "raw",
        "abelianization"=>abelianization_invariants(length(get(small,name,object).generator_matrices),get(small,name,object).relators)) for (name,object) in [(:GL,C),(:SL,S),(:PSL,Q)]))
    run["total_in_process_seconds"]=(time_ns()-start)/1e9;run["source_hash"]=ctx.source_hash;save()
    initial=joinpath(directory,"initial_complete_run.json")
    isfile(initial) || atomic_json(initial,run)
    println("D$D completed: ",run["geometry"],"; seconds ",run["total_in_process_seconds"]);flush(stdout)
    ArithmeticGroupResult(E,ctx,C,S,Q,small,metadata,abspath(directory),run)
end
