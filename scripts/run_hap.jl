# julia --project=. scripts/run_hap.jl --requests=examples/voronoi_requests.json
# Computes fresh complexes through HAP's general record interface, without Magma.
using ArithmeticVoronoi, Oscar, JSON3
import Oscar: GAP
include("request_arguments.jl")

function main(args)
    degrees=[parse(Int,split(arg,'=';limit=2)[2]) for arg in args if startswith(arg,"--degree=")]
    length(degrees)<=1 || throw(ArgumentError("specify --degree only once"))
    degree=isempty(degrees) ? 2 : only(degrees)
    degree>0 || throw(ArgumentError("homology degree must be positive"))
    (;root,resume,requests)=request_arguments(filter(arg->!startswith(arg,"--degree="),args);
        default_discriminants=[-7,5],default_output=joinpath(ArithmeticVoronoi.PACKAGE_ROOT,"results/hap_fields"))
    for r in requests
        r["group"]=="SL" && r["field"]["discriminant"]>0 &&
            throw(ArgumentError("real-quadratic SL complexes are not implemented"))
    end
    enable_hap_julia_backend(;progress=true,resolution_length=degree+1)
    # GAP.jl's callable wrapper keeps record options inside GAP.
    construct=GAP.evalstr("request -> VoronoiComplex(request : backend := \"julia\")")
    mkpath(root)
    entries=Dict{String,Any}()
    summary=Dict("schema"=>"av_hap_batch_v1","requests"=>requests,"degree"=>degree,
        "backend"=>"julia","status"=>"running","groups"=>entries,
        "julia_version"=>string(VERSION),"oscar_version"=>string(pkgversion(Oscar)),
        "hap_version"=>String(GAP.gap_to_julia(GAP.Globals.GAPInfo.PackagesLoaded.hap[2])))
    save()=ArithmeticVoronoi.atomic_json(joinpath(root,"summary.json"),summary)
    for request in requests
        name="$(request["group"])_D$(request["field"]["discriminant"])"
        path=joinpath(root,name*".json")
        if resume && isfile(path)
            previous=JSON3.read(read(path,String),Dict{String,Any})
            if get(previous,"status",nothing)=="complete" && matching_request(previous,request) &&
                get(previous,"degree",nothing)==degree
                entries[name]=previous;save();continue
            end
        end
        entry=Dict{String,Any}("request"=>request,"degree"=>degree,"status"=>"running")
        entries[name]=entry;save()
        try
            started=time_ns()
            K=construct(GAP.GapObj(request;recursive=true))
            entry["complex_seconds"]=(time_ns()-started)/1e9
            started=time_ns()
            R=GAP.Globals.FreeGResolution(K,degree+1)
            C=GAP.Globals.TensorWithIntegers(R)
            entry["homology"]=Dict(string(n)=>sort!(BigInt.(GAP.gap_to_julia(GAP.Globals.Homology(C,n)))) for n in 1:degree)
            entry["homology_seconds"]=(time_ns()-started)/1e9
            entry["status"]="complete"
        catch err
            err isa InterruptException && rethrow()
            entry["status"]="failed";entry["error"]=sprint(showerror,err)
            showerror(stderr,err,catch_backtrace());println(stderr)
        end
        ArithmeticVoronoi.atomic_json(path,entry);save()
    end
    summary["status"]=all(e->e["status"]=="complete",values(entries)) ? "complete" : "incomplete"
    save()
    summary["status"]=="complete" || error("some homology computations failed; see $(joinpath(root,"summary.json"))")
    println("Homology results written to ",abspath(root))
end

if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS)
end
