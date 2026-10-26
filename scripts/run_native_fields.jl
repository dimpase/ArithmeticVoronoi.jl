# Run from the checkout root with julia --project=. scripts/run_native_fields.jl.
# --output=DIR selects the parent directory; --resume retains completed fields.
# --requests=FILE reads JSON records. Integers are quadratic GL shortcuts.
const import_started=time_ns()
using ArithmeticVoronoi, Oscar, LinearAlgebra, JSON3
const NAV=ArithmeticVoronoi
const import_seconds=(time_ns()-import_started)/1e9
const SURVEY_DISCRIMINANTS=[5,8,12,13,17,21,24,28,29,33,37,40,41,44,56,60]
include("request_arguments.jl")

function validate_words(R;resume=false)
    P=R.enumeration.problem;K=P.field;w=P.generator
    e=R.sl.units.rank==0 ? R.sl.units.torsion_generator : R.sl.units.epsilon
    n=P.rank;I=identity_matrix(K,n)
    elementary(i,j,a)=begin g=deepcopy(I);g[i,j]=a;g end
    U=elementary(1,n,1);W=elementary(1,n,w)
    L=elementary(n,1,1);Lw=elementary(n,1,w)
    H=deepcopy(I);H[1,1]=e;H[n,n]=inv(e)
    J=deepcopy(I);J[1,1]=e
    inputs=[("minus_identity",-I),("U",U),("W",W),("lower_1",L),
        ("lower_w",Lw),("unit_diagonal_SL",H),("product",U*inv(Lw)*H*W),
        ("scalar",e*I),("inverse_scalar",inv(e)*I),
        ("unit_diagonal_GL",J)]
    R.sl.center_order>2 && push!(inputs,("SL_center",R.sl.center_generator*I))
    examples=Dict{String,Any}[]
    for (label,g) in inputs,target in (:GL,:SL,:PSL)
        target!=:GL && det(g)!=1 && continue
        path=joinpath(R.directory,"words","$(label)_$target.json")
        c=resume ? matrix_to_word(R,g;target,cache_path=path) : matrix_to_word(R,g;target)
        resume || NAV.atomic_json(path,c)
        NAV.ensure(verify_word_certificate(R,g,c),"word certificate failed for $label/$target")
        push!(examples,Dict("input"=>label,"group"=>string(target),
            "word_length"=>length(c["word"]),"certificate"=>relpath(path,R.directory)))
    end
    # Also exercise the original presentation alphabets independently of Tietze maps.
    for target in (:GL,:SL,:PSL)
        c=matrix_to_word(R,W;target,simplified=false)
        NAV.ensure(verify_word_certificate(R,W,c),"raw word certificate failed for $target")
        NAV.atomic_json(joinpath(R.directory,"words","W_$(target)_raw.json"),c)
    end
    validation=Dict("status"=>"passed","distinct_inputs"=>length(inputs),
        "simplified_certificates"=>length(examples),"raw_certificates"=>3,"examples"=>examples)
    NAV.atomic_json(joinpath(R.directory,"word_validation.json"),validation)
    R.run["status"]["constructive_words"]="validated_GL_SL_PSL"
    NAV.atomic_json(joinpath(R.directory,"run.json"),R.run)
    validation
end

field_key(request)=request["rank"]==2 ? string(request["field"]["discriminant"]) :
    "rank$(request["rank"])_D$(request["field"]["discriminant"])"
field_directory(request)=request["rank"]==2 ? "D$(request["field"]["discriminant"])" :
    "rank$(request["rank"])/D$(request["field"]["discriminant"])"

function compute_field(request,root;on_result=nothing,resume=false)
    D=request["field"]["discriminant"]
    path=field_directory(request)
    R=run_field_pipeline(request;backend=:native,directory=joinpath(root,path),
        resume)
    println("D$D: checking matrix-to-word certificates");flush(stdout)
    started=time_ns();validation=validate_words(R;resume);seconds=(time_ns()-started)/1e9
    comparison=isnothing(on_result) ? nothing : on_result(R)
    Dict("status"=>"complete","directory"=>path,"request"=>request,
        "perfect_form_orbits"=>length(R.enumeration.representatives),
        "labelled_neighbours"=>length(R.enumeration.neighbours),
        "oriented_loops"=>R.run["geometry"]["oriented_loops"],
        "ordinary_class_number"=>R.metadata.ordinary_class_number,
        "fundamental_unit"=>R.metadata.fundamental_unit,
        "fundamental_unit_norm"=>R.metadata.fundamental_unit_norm,
        "pipeline_seconds"=>R.run["total_in_process_seconds"],
        "stage_seconds"=>R.run["stage_seconds"],"presentations"=>R.run["presentations"],
        "word_validation_seconds"=>seconds,"word_validation"=>validation,"comparison"=>comparison)
end

function main(args;on_result=nothing,default_discriminants=SURVEY_DISCRIMINANTS,
              default_output=joinpath(pwd(),"results/native_fields"),validate_request=identity)
    (;root,resume,requests)=request_arguments(args;default_discriminants,default_output)
    requests=NAV.enumeration_request.(requests)
    foreach(validate_request,requests)
    discriminants=[r["field"]["discriminant"] for r in requests]
    LinearAlgebra.BLAS.set_num_threads(1)
    mkpath(root)
    fields=Dict{String,Any}()
    summary_path=joinpath(root,"summary.json")
    if resume && isfile(summary_path)
        previous=JSON3.read(read(summary_path,String),Dict{String,Any})
        NAV.ensure(previous["schema"]=="av_native_field_batch_v1","incompatible batch summary")
        for request in requests
            key=field_key(request)
            haskey(previous["fields"],key) || continue
            fields[key]=previous["fields"][key]
        end
    end
    summary=Dict("schema"=>"av_native_field_batch_v1","status"=>"running",
        "discriminants"=>discriminants,"requests"=>requests,"fields"=>fields,"resumed"=>resume,
        "julia_version"=>string(VERSION),"oscar_version"=>string(pkgversion(Oscar)),
        "julia_threads"=>Threads.nthreads(),"blas_threads"=>LinearAlgebra.BLAS.get_num_threads(),
        "package_import_seconds"=>import_seconds,
        "timing_note"=>"Fields run sequentially. The first field computed in each Julia session includes compilation; saved certificates are replayed only when --resume is supplied.")
    save()=NAV.atomic_json(summary_path,summary)
    save()
    for request in requests
        D=request["field"]["discriminant"];key=field_key(request)
        if resume && haskey(fields,key) && fields[key]["status"]=="complete" &&
            matching_request(fields[key],request) &&
            (isnothing(on_result) || !isnothing(get(fields[key],"comparison",nothing)))
            println(field_directory(request),": retaining completed batch result");flush(stdout)
            continue
        end
        fields[key]=Dict("status"=>"running","request"=>request);save()
        try
            fields[key]=compute_field(request,root;on_result,resume)
        catch err
            err isa InterruptException && rethrow()
            fields[key]=Dict("status"=>"failed","request"=>request,"error"=>sprint(showerror,err))
            showerror(stderr,err,catch_backtrace());println(stderr);flush(stderr)
        end
        save();GC.gc()
    end
    failed=[field_key(request) for request in requests if fields[field_key(request)]["status"]!="complete"]
    summary["status"]=isempty(failed) ? "complete" : "incomplete";save()
    isempty(failed) || error("Failed fields: $(join(failed,", ")); inspect each D*/run.json and the batch log")
    println("All $(length(discriminants)) fields completed; results in ",abspath(root))
    (;root,summary)
end

if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end
