module ArithmeticVoronoiData

using JSON3, SHA
export available_presentations, load_presentation

# This package contains recorded results only. It neither imports nor calls
# ArithmeticVoronoi, OSCAR, HAP or Magma. No data is read at package import.
const COLLECTION_FILES = (
    native_real="native_presentations.json",
    native_imaginary="imaginary_presentations.json",
    page="page_presentations.json",
)

function read_collection(collection::Symbol)
    haskey(COLLECTION_FILES,collection) || throw(ArgumentError(
        "unknown recorded collection; choose :native_real, :native_imaginary or :page"))
    path=joinpath(@__DIR__,"../data",COLLECTION_FILES[collection])
    raw=read(path)
    (;bundle=JSON3.read(String(copy(raw)),Dict{String,Any}),
      path=abspath(path),sha256=bytes2hex(sha256(raw)))
end

function entries(collection,bundle)
    collection==:page && return [
        (;discriminant=Int(g["discriminant"]),group=Symbol(g["group"]),data=g,field=nothing)
        for g in bundle["groups"]]
    [(;discriminant=parse(Int,D),group=Symbol(name),data=g,field=f["field"])
        for (D,f) in bundle["fields"] for (name,g) in f["groups"]]
end

"""List the recorded presentations in an explicitly selected collection."""
function available_presentations(collection::Symbol)
    source=read_collection(collection)
    sort!([(;collection,e.discriminant,rank=2,e.group)
        for e in entries(collection,source.bundle)];by=e->(e.discriminant,string(e.group)))
end

"""
    load_presentation(collection, discriminant; group=:GL, rank=2)

Load a PRECOMPUTED presentation from :native_real, :native_imaginary or :page.
This never runs an arithmetic computation. Missing entries raise ArgumentError.
The return value identifies the recorded source and its SHA-256 checksum.
Matrices use exact rational coefficient strings in the recorded quadratic basis.
"""
function load_presentation(collection::Symbol,D::Integer;group::Symbol=:GL,rank::Integer=2)
    rank==2 || throw(ArgumentError("this recorded-results library currently contains rank-two presentations only"))
    group in (:GL,:SL,:PSL) || throw(ArgumentError("choose group=:GL, :SL or :PSL"))
    source=read_collection(collection)
    matches=filter(e->e.discriminant==D && e.group==group,entries(collection,source.bundle))
    length(matches)==1 || throw(ArgumentError("no recorded $group presentation for discriminant $D in $collection"))
    entry=only(matches);data=entry.data
    matrices=collection==:native_real ? [g["matrix"] for g in data["generators"]] : data["generator_matrices"]
    abelianization=if collection==:page
        factors=data["abelian_invariants"]
        Dict("free_rank"=>count(iszero,factors),"torsion"=>string.(filter(x->x>1,factors)))
    else
        data["abelianization"]
    end
    @info "Loading a recorded presentation; no computation is performed" collection discriminant=D group source=source.path
    (;origin=:recorded_library,collection,discriminant=Int(D),rank=2,group,
      field=entry.field,generator_matrices=matrices,relators=[Int.(w) for w in data["relators"]],
      abelianization,projective=group==:PSL,
      provenance=(;package="ArithmeticVoronoiData",version=string(pkgversion(@__MODULE__)),
                  source_file=source.path,source_sha256=source.sha256,schema=source.bundle["schema"]))
end

end
