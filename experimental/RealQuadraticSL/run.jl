# julia --project=. experimental/RealQuadraticSL/run.jl 5 /durable/output/directory [rank]
include("RealQuadraticSL.jl")
using .RealQuadraticSL
using Oscar
import Oscar:GAP
D=isempty(ARGS) ? 5 : parse(Int,ARGS[1])
output=length(ARGS)>=2 ? abspath(ARGS[2]) : nothing
n=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 2
result=compute_real_sl_complex(D;rank=n,progress=true,certificate_directory=output)
R=GAP.Globals.FreeGResolution(result.hap,4)
C=GAP.Globals.TensorWithIntegers(R)
println("Integral SL",n," homology H₁, H₂, H₃: ",
    [GAP.gap_to_julia(GAP.Globals.Homology(C,n)) for n in 1:3])
println("Construction times in seconds: ",result.timings)
