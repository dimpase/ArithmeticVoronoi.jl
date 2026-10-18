using ArithmeticVoronoi, Oscar, JSON3
const ROOT=normpath(joinpath(@__DIR__,".."))
for (name,P) in [("classical2",classical_problem()),("D5",arithmetic_problem(5)),("D8",arithmetic_problem(8))]
    data=JSON3.read(read(joinpath(ROOT,"results/export_probe",name,"output.json"),String),Dict{String,Any})
    elapsed=@elapsed enumeration=validate_export(P,data)
    println(name,": ",enumeration.completion_status," reps=",length(enumeration.representatives)," neighbours=",length(enumeration.neighbours)," validation_seconds=",elapsed)
    println(enumeration.provenance["closure"])
end
