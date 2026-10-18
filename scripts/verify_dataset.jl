using ArithmeticVoronoi
directory=isempty(ARGS) ? joinpath(@__DIR__,"../results/dataset/v1") : ARGS[1]
level="--complete" in ARGS ? :complete : :algebraic
report=verify_dataset(directory;level)
println(report)
