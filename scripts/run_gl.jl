# julia --project=. scripts/run_gl.jl --requests=examples/higher_rank_gl_requests.json
using ArithmeticVoronoi, JSON3
include("request_arguments.jl")

function run_gl_batch(args)
    (;root,resume,requests)=request_arguments(args;default_discriminants=[-4],
        default_output=joinpath(ArithmeticVoronoi.PACKAGE_ROOT,"results/gl_presentations"))
    all(r->r["group"]=="GL",requests) || throw(ArgumentError("run_gl.jl requires GL requests"))
    # Enumeration is fresh. --resume explicitly enables exact replay of saved
    # attaching-loop certificates in the selected output directory.
    entries=Dict{String,Any}()
    for request in requests
        label="GL$(request["rank"])_D$(request["field"]["discriminant"])"
        R=run_gl_pipeline(request;directory=joinpath(root,label),resume)
        entries[label]=R.run
        ArithmeticVoronoi.atomic_json(joinpath(root,"summary.json"),entries)
    end
    println("GL presentations written to ",abspath(root))
end

if abspath(PROGRAM_FILE)==@__FILE__
    run_gl_batch(ARGS)
end
