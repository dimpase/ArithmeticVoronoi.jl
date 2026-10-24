# Run from the package root: julia --startup-file=no --project=. scripts/word_example.jl 12
using ArithmeticVoronoi, Oscar

D = isempty(ARGS) ? 12 : parse(Int, only(ARGS))
directory = joinpath(@__DIR__, "../results/word_example", "D$D")
R = run_field_pipeline(voronoi_request(D); backend=:polyhedral_common, directory)
P = R.enumeration.problem
g = matrix(P.field, [1 P.generator; 0 1])
path = joinpath(directory, "words", "tau_w.json")
c = matrix_to_word(R, g; target=:SL, cache_path=path)
@assert verify_word_certificate(R, g, c)
replayed = matrix_to_word(R, g; target=:SL, cache_path=path)
@assert replayed["word"] == c["word"]
println("Input matrix: ", g)
println("Word in simplified SL generators: ", c["word"])
println("Generator matrices: ", R.simplified[:SL].generator_matrices)
println("Verified certificate and cache replay: ", abspath(path))
