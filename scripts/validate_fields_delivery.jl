using ArithmeticVoronoi,Oscar,JSON3
include(joinpath(@__DIR__,"../test/runtests.jl"))
include(joinpath(@__DIR__,"validate_sanitizer.jl"))
directory=export_dataset(FIELD_WORD_RESULTS)
verification=verify_dataset(directory;level=:complete)
ArithmeticVoronoi.atomic_json(joinpath(@__DIR__,"../results/fields_dataset_complete_verification.json"),verification)
