using ArithmeticVoronoi, Oscar
const AV=ArithmeticVoronoi
P=arithmetic_problem(5)
println("Loading the validated D5 enumeration");flush(stdout)
E=enumerate_perfect_forms(P)
seconds=@elapsed context=TwoCellContext(E)
println("Exact ridges per cone: ",length.(context.ridges),"; context seconds ",seconds);flush(stdout)
seconds=@elapsed certificate=walk_ridge(context,1,1)
write_first_two_cell(context,certificate)
println("First two-cell: ",length(certificate.steps)," crossings; residual ",certificate.residual,
        "; residual word ",certificate.residual_word,"; seconds ",seconds)
