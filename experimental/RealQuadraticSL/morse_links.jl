using ArithmeticVoronoi, Oscar, JSON
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
    exponent_cache::Dict{Any,Int}
    flags::Any
    orbit_flags::Any
end

function morse_experiment(C;progress=false)
    P=C.enumeration.problem
    @assert P.degree==2 && P.discriminant>0 && P.rank==2
    units=AV.verified_units(P)
    @assert all(layer->all(S->all(g->iszero(AV.unit_coordinates(P,units,det(g)).exponent),S),layer),C.stabilizers) "height is not invariant under a cell stabilizer"
    I=identity_matrix(P.field,P.rank)
    vertices=[[MorseVertex(k-1,i,I,(k-1,AV.cell_vector_key(c.vectors)),0)
        for (i,c) in enumerate(layer)] for (k,layer) in enumerate(C.cells)]
    M=MorseExperiment(C,units,vertices,Dict(),Dict(one(P.field)=>0),
        [[Vector{MorseVertex}[] for c in layer] for layer in C.cells],
        [[Vector{MorseVertex}[] for c in layer] for layer in C.cells])
    for (k,layer) in enumerate(C.cells),(i,c) in enumerate(layer)
        top=vertices[k][i]
        if k==1
            M.flags[k][i]=[[top]]
        else
            for (j,g) in C.boundaries[k][i],flag in M.flags[k-1][abs(j)]
                push!(M.flags[k][i],[MorseVertex[morse_act(M,g,v) for v in flag];top])
            end
        end
        seen=Set{Any}()
        for flag in M.flags[k][i]
            key=Tuple(v.key for v in flag)
            key in seen && continue
            push!(M.orbit_flags[k][i],flag)
            for s in C.stabilizers[k][i]
                push!(seen,Tuple(morse_act(M,s,v).key for v in flag))
            end
        end
        if progress
            println("Flag orbits for cell ",(k-1,i),": ",length(M.orbit_flags[k][i])," / ",length(M.flags[k][i]));flush(stdout)
        end
    end
    M
end

function morse_act(M,g,v)
    P=M.complex.enumeration.problem
    key=(AV.matrix_key(AV.encode_matrix(P,g)),v.key)
    get!(M.action_cache,key) do
        t=g*v.translate
        d=det(t)
        shift=get!(M.exponent_cache,d) do
            Int(AV.unit_coordinates(P,M.units,d).exponent)
        end
        V=AV.cell_image(P,t,M.complex.cells[v.dim+1][v.orbit].vectors)
        MorseVertex(v.dim,v.orbit,t,(v.dim,AV.cell_vector_key(V)),shift)
    end
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
    w1=(P.root_sum+sqrt(Float64(P.discriminant)))/2
    w2=(P.root_sum-sqrt(Float64(P.discriminant)))/2
    at(x,w)=Float64(coeff(x,0))+w*Float64(coeff(x,1))
    epsilon1=at(M.units.epsilon,w1)
    bases=Float64[]
    for (k,layer) in enumerate(C.cells),(i,c) in enumerate(layer)
        F=sum(AV.act_form(P,c.point,inv(s)) for s in C.stabilizers[k][i])/length(C.stabilizers[k][i])
        delta=det(F)
        push!(bases,-log(at(delta,w1)/at(delta,w2))/(4log(epsilon1)))
    end
    data=Dict("D"=>P.discriminant,"base_heights"=>bases,"orbit_dimensions"=>[k-1 for (k,l) in enumerate(C.cells) for c in l],
        "links"=>[Dict("vertices"=>[Dict("orbit"=>ids[v.dim+1][v.orbit],"shift"=>v.shift) for v in link.vertices],
                       "facets"=>[collect(f) for f in link.facets]) for link in links])
    data,links
end

function simplicial_faces(facets)
    faces=Set{Tuple}()
    for facet in facets
        for mask in 1:(1<<length(facet))-1
            push!(faces,Tuple(facet[j] for j in eachindex(facet) if !iszero(mask & (1<<(j-1)))))
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
    if isnothing(bases)
        bases=[(round(Int,x*1000)*1_000_000+i)//1_000_000_000 for (i,x) in enumerate(data["base_heights"])]
    end
    @assert length(unique(mod.(bases,1)))==length(bases)
    reports=Any[]
    for (i,link) in enumerate(data["links"]),sgn in (-1,1)
        delta=[bases[v["orbit"]]+v["shift"]-bases[i] for v in link["vertices"]]
        @assert all(!iszero,delta)
        facets=[Tuple(v for v in f if sgn*delta[v]>0) for f in link["facets"]]
        filter!(!isempty,facets)
        result=collapse_link(facets)
        push!(reports,(;orbit=i,sgn,result))
        if verbose
            println("Link ",i," ",sgn,"; ",length(unique(v for f in facets for v in f))," vertices; collapsible=",result.passed,"; residual ",[count(f->length(f)==k,result.residual) for k in 1:4]);flush(stdout)
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
