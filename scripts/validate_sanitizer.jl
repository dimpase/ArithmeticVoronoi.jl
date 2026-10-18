using ArithmeticVoronoi, Oscar, JSON3
const ROOT=normpath(joinpath(@__DIR__,".."))
const AV=ArithmeticVoronoi
extended="--fields" in ARGS
relative=extended ? "results/sanitizer_fields" : "results/sanitizer"
path=joinpath(ROOT,relative,"summary.json")
summary=JSON3.read(read(path,String),Dict{String,Any})
@assert summary["status"]=="passed"
examples=extended ? [("D$D",arithmetic_problem(D)) for D in (12,13,40)] : [("classical2",classical_problem()),("D5",arithmetic_problem(5)),("D8",arithmetic_problem(8))]
for (name,P) in examples
    output=joinpath(ROOT,relative,name,"output.json")
    seconds=@elapsed checked=validate_export(P,JSON3.read(read(output,String),Dict{String,Any}))
    release=enumerate_perfect_forms(P)
    matches=[findall(r->!isnothing(equivalence(P,r.form,s.form)),release.representatives) for s in checked.representatives]
    @assert all(length(m)==1 for m in matches)
    @assert sort(only.(matches))==collect(eachindex(release.representatives))
    @assert length(checked.neighbours)==length(release.neighbours)
    summary["cases"][name]["math_validation"]="passed all exact invariants and matched release orbits K-linearly"
    summary["cases"][name]["math_validation_seconds"]=seconds
    summary["cases"][name]["orbit_map_to_release"]=only.(matches)
    println(name,": sanitizer output independently validated and equivalent to release")
end
summary["independent_math_validation"]="passed"
AV.atomic_json(path,summary)
