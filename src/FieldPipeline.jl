export ArithmeticGroupResult, run_field_pipeline, load_simplification

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

function load_simplification(P,path,matrices,relations)
    d=JSON3.read(read(path,String),Dict{String,Any})
    S=(;label=d["group"],projective=get(d,"projective",false),algorithm=d["algorithm"],
      isomorphism_verified=d["isomorphism_verified"],old_generator_count=d["old_generator_count"],
      old_relator_count=d["old_relator_count"],gap_version=d["gap_version"],
      generator_matrices=Any[field_matrix_wire(P,g["matrix"]) for g in d["generators"]],
      backward_generator_words=[Int.(g["word_in_old_generators"]) for g in d["generators"]],
      forward_generator_words=[Int.(w) for w in d["old_generators_as_words_in_new"]],
      relators=[Int.(w) for w in d["relators"]])
    verify_simplification(matrices,relations,S)
    S
end

function run_field_pipeline(D;directory=joinpath(PACKAGE_ROOT,"results/fields","D$D"),
                            simplify=true,timeout=600.0,reuse_baseline=D in (5,8))
    mkpath(directory);start=time_ns();timings=Dict{String,Float64}()
    if isfile(joinpath(directory,"run.json"))
        history=joinpath(directory,"history");mkpath(history)
        cp(joinpath(directory,"run.json"),joinpath(history,"run_$(time_ns()).json"))
        previous=JSON3.read(read(joinpath(directory,"run.json"),String),Dict{String,Any})
        initial=joinpath(directory,"initial_complete_run.json")
        if !isfile(initial) && haskey(previous,"total_in_process_seconds") && all(get(previous["status"],s,"pending")=="complete" for s in ("GL","SL","PSL"))
            cp(joinpath(directory,"run.json"),initial)
        end
    end
    status=Dict{String,Any}("enumeration"=>"pending","two_cells"=>"pending","GL"=>"pending",
        "SL"=>"pending","PSL"=>"pending","simplification"=>Dict(k=>(simplify ? "pending" : "not_requested") for k in ("GL","SL","PSL")),"constructive_words"=>"unvalidated",
        "literature_comparison"=>D==40 ? "not_applicable" : "pending")
    run=Dict{String,Any}("schema"=>"av_field_run_v1","discriminant"=>D,"status"=>status,
        "stage_seconds"=>timings,"julia_version"=>string(VERSION),"oscar_version"=>string(pkgversion(Oscar)),
        "julia_threads"=>Threads.nthreads(),"started_utc"=>string(now(UTC)))
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
    P=arithmetic_problem(D)
    metadata=stage("field_and_unit_verification",()->field_metadata(P));atomic_json(joinpath(directory,"field.json"),metadata)
    E=stage("perfect_enumeration_and_validation",()->enumerate_perfect_forms(P;timeout))
    status["enumeration"]="complete";run["enumeration_from_cache"]=E.provenance["from_cache"]
    run["backend_seconds"]=get(E.provenance,"backend_seconds",nothing)
    write_enumeration(joinpath(directory,"enumeration.json"),E);save()
    ctx=stage("ridges",()->TwoCellContext(E))
    certdir=reuse_baseline ? joinpath(PACKAGE_ROOT,"results/two_cells","D$D") : joinpath(directory,"two_cells")
    cells=stage("two_cells_and_coverage",()->begin
        ledgerpath=joinpath(certdir,"coverage.json")
        saved=isfile(ledgerpath) ? JSON3.read(read(ledgerpath,String),Dict{String,Any}) : nothing
        if reuse_baseline || (!isnothing(saved) && saved["status"]=="complete")
            ledger=saved
            ensure(ledger["source_hash"]==ctx.source_hash,"complete ridge cache belongs to another enumeration")
            certificates=TwoCellCertificate[read_two_cell(ctx,joinpath(certdir,file)) for e in ledger["entries"] for file in e["certificates"]]
            for k in 1:2:length(certificates);reverse_cycle_check(ctx,certificates[k],certificates[k+1]);end
            (;certificates,ledger)
        else
            enumerate_two_cells(ctx;directory=certdir)
        end
    end)
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
        original=joinpath(PACKAGE_ROOT,"results/two_cells","D$D",string(name)*"_simplified.json")
        small[name]=stage(string(name)*"_simplification",()->begin
            if reuse_baseline && name!=:PSL
                load_simplification(P,original,object.generator_matrices,object.relators)
            elseif isfile(path)
                load_simplification(P,path,object.generator_matrices,object.relators)
            else
                simplify_matrix_presentation(object.generator_matrices,object.relators;label="$name(2,O), D$D",projective=name==:PSL)
            end
        end)
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
