include(joinpath(@__DIR__,"../integrations/perfect_forms/adapter.jl"))
using .PerfectFormsAdapter
using ArithmeticVoronoi
for D in (isempty(ARGS) ? [12,13,40] : parse.(Int,ARGS))
    run_field_pipeline(voronoi_request(D);backend=:polyhedral_common)
end
