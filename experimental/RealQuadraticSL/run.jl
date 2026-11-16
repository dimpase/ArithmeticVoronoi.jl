# julia --project=. experimental/RealQuadraticSL/run.jl 5 /durable/output/directory
include("RealQuadraticSL.jl")
using .RealQuadraticSL
using Oscar
import Oscar:GAP
D=isempty(ARGS) ? 5 : parse(Int,ARGS[1])
output=length(ARGS)>=2 ? abspath(ARGS[2]) : nothing
result=compute_real_sl_complex(D;progress=true,certificate_directory=output)
R=GAP.Globals.FreeGResolution(result.hap,4)
C=GAP.Globals.TensorWithIntegers(R)
println("Integral SL₂ homology H₁, H₂, H₃: ",
    [GAP.gap_to_julia(GAP.Globals.Homology(C,n)) for n in 1:3])
println("Construction times in seconds: ",result.timings)
