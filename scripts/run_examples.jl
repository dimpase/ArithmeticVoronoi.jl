include(joinpath(@__DIR__,"../integrations/perfect_forms/adapter.jl"))
using .PerfectFormsAdapter
using ArithmeticVoronoi, Oscar, JSON3
import Hecke
const ROOT=normpath(joinpath(@__DIR__,".."))
const AV=ArithmeticVoronoi
mkpath(joinpath(ROOT,"results"))
backend=default_backend()
report=Dict{String,Any}("julia_version"=>string(VERSION),"oscar_version"=>string(pkgversion(Oscar)),
    "hecke_version"=>string(pkgversion(Hecke)),"gap_jl_version"=>string(pkgversion(Oscar.GAP)),
    "gap_version"=>string(Oscar.GAP.gap_to_julia(Oscar.GAP.Globals.GAPInfo.Version)),
    "executable"=>backend.executable,"executable_sha256"=>backend.executable_sha256,
    "backend_commit"=>backend.commit,"examples"=>Any[])
for (name,P) in [("classical2",classical_problem()),("D5",arithmetic_problem(5)),("D8",arithmetic_problem(8))]
    println("START ",name);flush(stdout)
    cold_seconds=@elapsed cold=enumerate_perfect_forms(P;backend=:polyhedral_common,config=backend,force=true)
    write_enumeration(joinpath(ROOT,"results",name*".json"),cold)
    warm_seconds=@elapsed warm=enumerate_perfect_forms(P;backend=:polyhedral_common,config=backend,force=false)
    @assert warm.completion_status==cold.completion_status==:complete
    @assert [r.trace_gram for r in warm.representatives]==[r.trace_gram for r in cold.representatives]
    @assert warm.provenance["from_cache"] && !cold.provenance["from_cache"]
    push!(report["examples"],Dict("name"=>name,"status"=>string(cold.completion_status),
        "orbits"=>length(cold.representatives),"stabilizer_orders"=>string.([r.stabilizer_order for r in cold.representatives]),
        "minimal_vector_pair_counts"=>[length(r.minimal_vectors) for r in cold.representatives],
        "facets"=>[length(r.facets) for r in cold.representatives],
        "facet_orbits"=>cold.provenance["closure"],"labelled_neighbours"=>length(cold.neighbours),
        "cold_in_process_seconds"=>cold_seconds,"warm_in_process_seconds"=>warm_seconds,
        "cold_provenance"=>cold.provenance,"warm_provenance"=>warm.provenance))
    AV.atomic_json(joinpath(ROOT,"results/examples_summary.json"),report)
    println("DONE ",name,": ",length(cold.representatives)," orbits, ",length(cold.neighbours)," labelled facet orbits; cold/warm ",cold_seconds,"/",warm_seconds," seconds");flush(stdout)
end
