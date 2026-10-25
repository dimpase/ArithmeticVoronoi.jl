# Default suite: native computations and small, ordinary-Git test fixtures.
# No data repository, LFS payloads, HAP or C++ executable is required.
using Test

@testset "ArithmeticVoronoi native engine" begin
    include("engine_isolation.jl")
    include("checkpoint_opt_in.jl")
    include("model.jl")
    include("requests.jl")
    include("quadratic_units.jl")
    include("native.jl")
    include("imaginary.jl")
    include("gl_presentations.jl")
    test_gl_presentations()
end
