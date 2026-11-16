using ..ArithmeticVoronoi
using Oscar, JSON
const AV=ArithmeticVoronoi

struct MorseVertex
    dim::Int
    orbit::Int
    translate::Any
    key::Any
    shift::Int
end
mutable struct MorseExperiment
    complex::Any
    units::Any
    vertices::Any
    action_cache::Dict{Any,MorseVertex}
    vertex_cache::Dict{Any,MorseVertex}
    matrix_cache::IdDict{Any,Any}
    exponent_cache::Dict{Any,Int}
    full_flag_counts::Vector{Vector{BigInt}}
    orbit_flags::Any
end

# Keys are exact integers. Keep ordinary-sized coordinates inline; retaining a
# heap-allocated BigInt for each small coordinate dominates large face posets.
compact_integer(x)=typemin(Int)<=x<=typemax(Int) ? Int(x) : x
morse_vertex_key(dim,V)=(dim,Tuple(Tuple(compact_integer(x) for x in v) for v in AV.cell_vector_key(V)))

function morse_experiment(C;progress=false,maximal_flags=true)
    P=C.enumeration.problem
    P.degree==2 && P.discriminant>0 && C.group==:GL ||
        throw(ArgumentError("the determinant-height construction requires a real quadratic GL complex"))
    units=AV.verified_units(P)
    @assert all(layer->all(S->all(g->iszero(AV.unit_coordinates(P,units,det(g)).exponent),S),layer),C.stabilizers) "height is not invariant under a cell stabilizer"
    I=identity_matrix(P.field,P.rank)
    vertices=[[MorseVertex(k-1,i,I,morse_vertex_key(k-1,c.vectors),0)
        for (i,c) in enumerate(layer)] for (k,layer) in enumerate(C.cells)]
    interned=Dict{Any,MorseVertex}(v.key=>v for layer in vertices for v in layer)
    M=MorseExperiment(C,units,vertices,Dict(),interned,IdDict(),Dict(one(P.field)=>0),
        [zeros(BigInt,length(layer)) for layer in C.cells],
        [[Vector{MorseVertex}[] for c in layer] for layer in C.cells])
    for (k,layer) in enumerate(C.cells),(i,c) in enumerate(layer)
        M.full_flag_counts[k][i]=k==1 ? big(1) :
            sum(M.full_flag_counts[k-1][abs(j)] for (j,g) in C.boundaries[k][i];init=big(0))
        maximal_flags || continue
        flags,covered=morse_flag_orbits(M,vertices[k][i],C.stabilizers[k][i])
        @assert covered==M.full_flag_counts[k][i] "flag orbits do not cover the face lattice"
        M.orbit_flags[k][i]=flags
        if progress
            println("Flag orbits for cell ",(k-1,i),": ",length(flags)," / ",covered);flush(stdout)
        end
    end
    M
end

# Descend through one facet from each orbit under the stabilizer of the flag
# chosen so far. Restrict that stabilizer at each step. This enumerates maximal
# flags directly up to symmetry, without constructing all their translates.
function morse_flag_orbits(M,top,stabilizer)
    result=Vector{MorseVertex}[]
    covered=big(0)
    function descend(flag,S)
        v=first(flag)
        if v.dim==0
            push!(result,flag)
            @assert length(stabilizer)%length(S)==0
            covered+=div(length(stabilizer),length(S))
            return
        end
        facets=[morse_act(M,v.translate*g,M.vertices[v.dim][abs(j)])
                for (j,g) in M.complex.boundaries[v.dim+1][v.orbit]]
        keys=Set(w.key for w in facets)
        @assert length(keys)==length(facets) "a cell boundary repeats a facet"
        seen=Set{Any}()
        for w in facets
            w.key in seen && continue
            next=Any[]
            for s in S
                image=morse_act(M,s,w)
                @assert image.key in keys "flag stabilizer does not preserve its facets"
                push!(seen,image.key)
                image.key==w.key && push!(next,s)
            end
            descend([w;flag],next)
        end
    end
    descend([top],stabilizer)
    result,covered
end

function clear_action_cache!(M)
    empty!(M.action_cache)
    empty!(M.matrix_cache)
end

# Most incidence coordinates are small. Checked machine arithmetic has a
# BigInt fallback, including sign normalization of typemin(Int); no overflow
# or rounding can change a cell key.
function morse_ray_image(rows,small,v)
    if small && all(x->x isa Int,v)
        try
            w=Tuple(begin
                value=0
                for j in eachindex(v)
                    value=Base.Checked.checked_add(value,Base.Checked.checked_mul(row[j],v[j]))
                end
                value
            end for row in rows)
            i=findfirst(!iszero,w)
            @assert !isnothing(i) "a group matrix mapped a minimum vector to zero"
            return w[i]<0 ? map(Base.Checked.checked_neg,w) : w
        catch error
            error isa OverflowError || rethrow()
        end
    end
    w=Tuple(sum(BigInt(row[j])*BigInt(v[j]) for j in eachindex(v)) for row in rows)
    i=findfirst(!iszero,w)
    @assert !isnothing(i) "a group matrix mapped a minimum vector to zero"
    w[i]<0 && (w=map(-,w))
    map(compact_integer,w)
end

function morse_act(M,g,v)
    g==M.vertices[1][1].translate && return v # identity
    P=M.complex.enumeration.problem
    action=get!(M.matrix_cache,g) do
        A=AV.encode_matrix(P,g)
        rows=Tuple(Tuple(compact_integer(BigInt(ZZ(A[i,j]))) for j in 1:ncols(A)) for i in 1:nrows(A))
        (rows=rows,small=all(row->all(x->x isa Int,row),rows),
            key=Tuple(compact_integer(x) for x in AV.matrix_key(A)))
    end
    get!(M.action_cache,(action.key,v.key)) do
        # Act on the already-translated minimum set. Multiplying two field
        # matrices and recomputing a unit exponent is needed only for a new cell.
        image_key=(v.dim,Tuple(sort!([morse_ray_image(action.rows,action.small,w) for w in v.key[2]])))
        image=get!(M.vertex_cache,image_key) do
            d=det(g)
            shift=v.shift+get!(M.exponent_cache,d) do
                Int(AV.unit_coordinates(P,M.units,d).exponent)
            end
            MorseVertex(v.dim,v.orbit,g*v.translate,image_key,shift)
        end
        @assert image.orbit==v.orbit
        image
    end
end

# This numerical value only proposes a rational height. Exact link collapses
# certify the result. For a totally positive delta=A+B*sqrt(D), evaluate the
# larger conjugate without cancellation and recover the other via norm(delta).
function morse_base_height(P,units,F)
    bf(x)=BigFloat(BigInt(numerator(x)))/BigFloat(BigInt(denominator(x)))
    root=sqrt(BigFloat(P.discriminant))
    delta=det(F)
    a,b=coeff(delta,0),coeff(delta,1)
    middle=a+b*P.root_sum/2
    absolute_norm=norm(delta)
    @assert middle>0 && absolute_norm>0 "determinant is not totally positive"
    larger=bf(middle)+abs(bf(b))*root/2
    log_ratio=iszero(b) ? zero(BigFloat) : Int(sign(b))*(2log(larger)-log(bf(absolute_norm)))
    a,b=coeff(units.epsilon,0),coeff(units.epsilon,1)
    epsilon1=bf(a+b*P.root_sum/2)+bf(b)*root/2
    result=Float64(-log_ratio/(4log(epsilon1)))
    isfinite(result) || throw(ArgumentError("could not propose a finite determinant height"))
    result
end

function morse_link_data(M)
    C=M.complex;P=C.enumeration.problem
    ids=[[sum(length,C.cells[1:k-1];init=0)+i for i in eachindex(layer)] for (k,layer) in enumerate(C.cells)]
    links=[(vertices=MorseVertex[],indices=Dict{Any,Int}(),facets=Set{Tuple}()) for layer in C.cells for c in layer]
    for layer in M.orbit_flags,flags in layer,flag in flags,(j,v) in enumerate(flag)
        id=ids[v.dim+1][v.orbit]
        link=links[id]
        for s in C.stabilizers[v.dim+1][v.orbit]
            t=s*inv(v.translate)
            facet=Int[]
            for (l,w) in enumerate(flag)
                l==j && continue
                x=morse_act(M,t,w)
                k=get!(link.indices,x.key) do
                    push!(link.vertices,x);length(link.vertices)
                end
                push!(facet,k)
            end
            isempty(facet) || push!(link.facets,Tuple(sort(facet)))
        end
    end
    # The numerical determinant height only proposes a rational Morse function.
    # The later exact link-collapse checks, not these approximations, certify it.
    bases=Float64[]
    for (k,layer) in enumerate(C.cells),(i,c) in enumerate(layer)
        F=sum(AV.act_form(P,c.point,inv(s)) for s in C.stabilizers[k][i])/length(C.stabilizers[k][i])
        push!(bases,morse_base_height(P,M.units,F))
    end
    data=Dict("D"=>P.discriminant,"base_heights"=>bases,"orbit_dimensions"=>[k-1 for (k,l) in enumerate(C.cells) for c in l],
        "links"=>[Dict("vertices"=>[Dict("orbit"=>ids[v.dim+1][v.orbit],"shift"=>v.shift) for v in link.vertices],
                       "facets"=>[collect(f) for f in link.facets]) for link in links])
    data,links
end

# Avoid machine-word overflow when the flag length exceeds the Int bit width.
function subset_masks(n;include_empty=false)
    n>=0 || throw(ArgumentError("a subset size must be nonnegative"))
    stop=n<8sizeof(Int)-1 ? (1<<n)-1 : (big(1)<<n)-1
    (include_empty ? 0 : 1):stop
end

function simplicial_faces(facets)
    faces=Set{Tuple}()
    for facet in facets
        for mask in subset_masks(length(facet))
            push!(faces,Tuple(facet[j] for j in eachindex(facet) if isodd(mask>>(j-1))))
        end
    end
    faces
end
function collapse_link(facets)
    faces=simplicial_faces(facets)
    up=Dict(f=>Set{Tuple}() for f in faces)
    for f in faces
        length(f)==1 && continue
        for j in eachindex(f);push!(up[Tuple(f[k] for k in eachindex(f) if k!=j)],f);end
    end
    free=sort!(Tuple[f for f in faces if length(up[f])==1];by=f->(length(f),f))
    pairs=Tuple{Tuple,Tuple}[]
    while !isempty(free)
        f=pop!(free)
        f in faces && length(up[f])==1 || continue
        g=only(up[f]);@assert isempty(up[g])
        push!(pairs,(f,g))
        for cell in (g,f)
            delete!(faces,cell)
            if length(cell)>1
                for j in eachindex(cell)
                    h=Tuple(cell[k] for k in eachindex(cell) if k!=j)
                    delete!(up[h],cell)
                    length(up[h])==1 && push!(free,h)
                end
            end
        end
    end
    (;passed=length(faces)==1 && length(only(faces))==1,pairs,residual=faces)
end
function evaluate_morse_links(data;bases=nothing,verbose=true)
    bases=rational_morse_heights(data;bases)
    reports=Any[]
    for (i,link) in enumerate(data["links"]),sgn in (-1,1)
        delta=[bases[v["orbit"]]+v["shift"]-bases[i] for v in link["vertices"]]
        @assert all(!iszero,delta)
        facets=[Tuple(v for v in f if sgn*delta[v]>0) for f in link["facets"]]
        filter!(!isempty,facets)
        result=collapse_link(facets)
        push!(reports,(;orbit=i,sgn,result))
        if verbose
            println("Link ",i," ",sgn,"; ",length(unique(v for f in facets for v in f))," vertices; collapsible=",result.passed,"; residual ",[count(f->length(f)==k,result.residual) for k in 1:maximum(length,result.residual;init=0)]);flush(stdout)
        end
    end
    (;passed=all(x->x.result.passed,reports),bases,reports)
end

using Random
function randomized_collapse(facets,seed)
    faces=simplicial_faces(facets)
    up=Dict(f=>Set{Tuple}() for f in faces)
    for f in faces
        length(f)==1 && continue
        for j in eachindex(f);push!(up[Tuple(f[k] for k in eachindex(f) if k!=j)],f);end
    end
    rng=MersenneTwister(seed)
    pairs=Tuple{Tuple,Tuple}[]
    while true
        free=sort!(Tuple[f for f in faces if length(up[f])==1];by=f->(length(f),f))
        isempty(free) && break
        f=rand(rng,free);g=only(up[f]);@assert isempty(up[g])
        push!(pairs,(f,g))
        for cell in (g,f)
            delete!(faces,cell)
            if length(cell)>1
                for j in eachindex(cell)
                    h=Tuple(cell[k] for k in eachindex(cell) if k!=j)
                    delete!(up[h],cell)
                end
            end
        end
    end
    (;passed=length(faces)==1 && length(only(faces))==1,pairs,residual=faces)
end

# Independent replay uses the full face set, not the search's coface table.
function replay_collapse(facets,pairs)
    faces=simplicial_faces(facets)
    for (sigma,tau) in pairs
        @assert sigma in faces && tau in faces
        @assert length(tau)==length(sigma)+1 && issubset(sigma,tau)
        @assert Set(f for f in faces if issubset(sigma,f))==Set([sigma,tau])
        delete!(faces,sigma);delete!(faces,tau)
    end
    @assert length(faces)==1 && length(only(faces))==1
    true
end

function certify_morse_links(data;bases=nothing)
    initial=evaluate_morse_links(data;bases,verbose=false)
    reports=Any[]
    for r in initial.reports
        link=data["links"][r.orbit]
        delta=[initial.bases[v["orbit"]]+v["shift"]-initial.bases[r.orbit] for v in link["vertices"]]
        facets=[Tuple(v for v in f if r.sgn*delta[v]>0) for f in link["facets"]]
        filter!(!isempty,facets)
        result=r.result
        if !result.passed
            for seed in 1:32
                result=randomized_collapse(facets,seed)
                result.passed && break
            end
        end
        @assert result.passed "No collapse certificate for link $(r.orbit), direction $(r.sgn)"
        @assert replay_collapse(facets,result.pairs)
        push!(reports,(;orbit=r.orbit,sgn=r.sgn,result))
    end
    (;passed=true,bases=initial.bases,reports)
end
