using ArithmeticVoronoi
results=Dict(D=>run_field_pipeline(D) for D in (5,8,12,13,40))
include(joinpath(@__DIR__,"../test/fields_words.jl"))
run_fields_words_tests(results)
directory=export_dataset(results)
println(verify_dataset(directory;level=:complete))
