export voronoi_request, normalize_voronoi_request

const VoronoiRequestRecord = Union{NamedTuple,AbstractDict}

"""
    voronoi_request(D; group=:GL, rank=2)

Describe a group over the maximal order of the quadratic field with
fundamental discriminant `D`, on the free lattice O^rank. The returned record can
be passed to `enumerate_perfect_forms` (GL), `hap_voronoi_complex` (GL, or
imaginary-quadratic SL), and serialized as JSON. `run_field_pipeline` still
requires rank-two GL.
"""
function voronoi_request(D::Integer;group=:GL,rank=2)
    normalize_voronoi_request((;group,rank,
        field=(;kind="quadratic",discriminant=D),
        order=(;kind="maximal"),lattice=(;kind="free")))
end

function request_record(value, context, allowed, required=allowed)
    value isa VoronoiRequestRecord || throw(ArgumentError("$context must be a record"))
    all(k->k isa Union{Symbol,AbstractString},keys(value)) ||
        throw(ArgumentError("$context keys must be strings or symbols"))
    record=Dict{String,Any}(string(k)=>v for (k,v) in pairs(value))
    length(record)==length(value) || throw(ArgumentError("$context contains duplicate keys"))
    extra=setdiff(collect(keys(record)),allowed)
    isempty(extra) || throw(ArgumentError("unsupported $context fields: $(join(sort!(extra), ", "))"))
    all(k->haskey(record,k),required) || throw(ArgumentError("$context requires $(join(required, ", "))"))
    record
end

function request_integer(value, context)
    value isa Integer && !(value isa Bool) && typemin(Int)<=value<=typemax(Int) ||
        throw(ArgumentError("$context must be an integer representable as Int"))
    Int(value)
end

"""
    normalize_voronoi_request(request)

Validate a Julia/JSON request against the arithmetic currently implemented by
ArithmeticVoronoi. HAP's dispatcher accepts more general records; this provider
rejects unsupported fields, orders, lattices and ranks instead of ignoring them.
Quadratic fields accept exactly one of `discriminant` and `radicand`; the result
uses a fundamental discriminant. This does not construct a field or enumerate.
"""
function normalize_voronoi_request(request::VoronoiRequestRecord)
    r=request_record(request,"request",["group","rank","field","order","lattice"])
    r["group"] in ("GL","SL",:GL,:SL) || throw(ArgumentError("the Julia provider supports group GL or SL"))
    n=request_integer(r["rank"],"rank")
    n>=2 || throw(ArgumentError("the Julia Voronoi provider requires rank at least 2"))
    f=request_record(r["field"],"field",["kind","discriminant","radicand"],["kind"])
    f["kind"]=="quadratic" || throw(ArgumentError("the Julia provider currently requires a quadratic field"))
    xor(haskey(f,"discriminant"),haskey(f,"radicand")) ||
        throw(ArgumentError("quadratic field requires exactly one of discriminant or radicand"))
    D=if haskey(f,"radicand")
        d=request_integer(f["radicand"],"field.radicand")
        d!=0 && d!=1 && is_squarefree(ZZ(d)) ||
            throw(ArgumentError("field.radicand must define a squarefree quadratic extension"))
        request_integer(mod(d,4)==1 ? BigInt(d) : 4BigInt(d),"field.discriminant")
    else
        request_integer(f["discriminant"],"field.discriminant")
    end
    (D<0 || D>1) && Hecke.is_fundamental_discriminant(D) ||
        throw(ArgumentError("a fundamental quadratic discriminant is required"))
    o=request_record(r["order"],"order",["kind"])
    o["kind"]=="maximal" || throw(ArgumentError("the Julia provider currently requires the maximal order"))
    l=request_record(r["lattice"],"lattice",["kind"])
    l["kind"]=="free" || throw(ArgumentError("the Julia provider currently requires the free lattice O^rank"))
    Dict{String,Any}("group"=>string(r["group"]),"rank"=>n,
        "field"=>Dict{String,Any}("kind"=>"quadratic","discriminant"=>D),
        "order"=>Dict("kind"=>"maximal"),"lattice"=>Dict("kind"=>"free"))
end

function arithmetic_problem(request::VoronoiRequestRecord)
    r=normalize_voronoi_request(request)
    arithmetic_problem(r["field"]["discriminant"];rank=r["rank"])
end

function enumeration_request(request)
    r=normalize_voronoi_request(request)
    r["group"]=="GL" || throw(ArgumentError("perfect-form enumeration and the presentation pipeline require a GL request; the pipeline derives SL and PSL presentations too"))
    r
end

function enumerate_perfect_forms(request::VoronoiRequestRecord;backend=:native,kwargs...)
    r=enumeration_request(request)
    enumerate_perfect_forms(arithmetic_problem(r);backend,kwargs...)
end
