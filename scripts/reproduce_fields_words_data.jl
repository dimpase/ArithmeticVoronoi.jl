include(joinpath(@__DIR__,"../integrations/perfect_forms/adapter.jl"))
using .PerfectFormsAdapter
using ArithmeticVoronoi
results=Dict(D=>run_field_pipeline(voronoi_request(D);backend=:polyhedral_common) for D in (5,8,12,13,40))
include(joinpath(@__DIR__,"../test/fields_words.jl"))
run_fields_words_tests(results;report_path=joinpath(@__DIR__,"../results/fields_words_tests.json"))
directory=export_dataset(results)
println(verify_dataset(directory;level=:complete))
