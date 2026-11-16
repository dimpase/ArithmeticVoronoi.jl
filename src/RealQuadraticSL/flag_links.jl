# Links in an order complex are flag complexes. A vertex or edge whose link
# is a cone can be removed by elementary simplicial collapses. The cone test
# uses only the graph; it never materializes the exponentially many cliques.
# Boissonnat--Pritam, SoCG 2020, doi:10.4230/LIPIcs.SoCG.2020.19.

function morse_cell_faces(M;progress=false)
    C=M.complex
    last_report=time()
    faces=[[MorseVertex[] for c in layer] for layer in C.cells]
    for k in 2:length(C.cells),i in eachindex(C.cells[k])
        seen=Dict{Any,MorseVertex}()
        for (j,g) in C.boundaries[k][i]
            for v in [M.vertices[k-1][abs(j)];faces[k-1][abs(j)]]
                w=morse_act(M,g,v)
                seen[w.key]=w
            end
        end
        faces[k][i]=sort!(collect(values(seen));by=v->v.key)
        clear_action_cache!(M)
        if progress && time()-last_report>=10
            println("Cell face posets: dimension ",k-1,", orbit ",i,"/",length(C.cells[k]),
                "; ",length(seen)," proper faces");flush(stdout)
            last_report=time()
        end
    end
    faces
end

# Materialize one link graph at a time. Large levels have millions of faces;
# keeping every link's edge list and adjacency table simultaneously is wasteful.
struct MorseLinkGraphs <: AbstractVector{Dict{String,Any}}
    links::Vector{Dict{Any,MorseVertex}}
    ids::Vector{Vector{Int}}
    experiment::Union{Nothing,MorseExperiment}
end
MorseLinkGraphs(links,ids)=MorseLinkGraphs(links,ids,nothing)

function morse_link_symmetries(M,k,i,vertices)
    P=M.complex.enumeration.problem
    I=identity_matrix(P.field,P.rank)
    generators=AV.cell_generators(P,M.complex.stabilizers[k][i])
    filter!(g->g!=I && g!=-I,generators)
    isempty(generators) && return Vector{Int}[]
    positions=Dict(v.key=>j for (j,v) in enumerate(vertices))
    permutations=Vector{Int}[]
    for g in generators
        A=AV.encode_matrix(P,g)
        rows=Tuple(Tuple(compact_integer(BigInt(ZZ(A[a,b]))) for b in 1:ncols(A)) for a in 1:nrows(A))
        small=all(row->all(x->x isa Int,row),rows)
        rays=Dict{Any,Any}()
        # A ray occurs in many cells. Transform each distinct ray once.
        permutation=[positions[(v.dim,Tuple(sort!([get!(rays,ray) do
            morse_ray_image(rows,small,ray)
        end for ray in v.key[2]])))] for v in vertices]
        @assert all(vertices[permutation[j]].orbit==v.orbit &&
                    vertices[permutation[j]].shift==v.shift for (j,v) in enumerate(vertices))
        push!(permutations,permutation)
    end
    permutations
end

Base.size(graphs::MorseLinkGraphs)=(length(graphs.links),)
Base.IndexStyle(::Type{MorseLinkGraphs})=IndexLinear()
function Base.getindex(graphs::MorseLinkGraphs,i::Int)
    @boundscheck checkbounds(graphs.links,i)
    vertices=sort!(collect(values(graphs.links[i]));by=v->v.key)
    # Incidence is reverse inclusion of complete minimum sets. An inverted
    # ray index finds supersets without quadratic all-pairs comparisons.
    containing=Dict{Any,BitSet}()
    for (j,v) in enumerate(vertices),ray in v.key[2]
        list=get!(containing,ray) do;BitSet();end
        push!(list,j)
    end
    sizes=Dict(ray=>length(list) for (ray,list) in containing)
    edges=Tuple{Int,Int}[]
    for (b,v) in enumerate(vertices)
        rays=sort!(collect(v.key[2]);by=ray->sizes[ray])
        candidates=copy(containing[first(rays)])
        for ray in rays[2:end];intersect!(candidates,containing[ray]);end
        for a in candidates
            vertices[a].dim<v.dim && push!(edges,(a,b))
        end
    end
    link=Dict{String,Any}("vertices"=>[Dict("orbit"=>graphs.ids[v.dim+1][v.orbit],"shift"=>v.shift) for v in vertices],"edges"=>edges)
    if !isnothing(graphs.experiment)
        k=findfirst(layer->!isempty(layer) && first(layer)<=i<=last(layer),graphs.ids)
        index=i-first(graphs.ids[k])+1
        link["colours"]=[v.dim for v in vertices]
        link["symmetry_generators"]=morse_link_symmetries(graphs.experiment,k,index,vertices)
    end
    link
end

function morse_link_graphs(M;progress=false)
    C=M.complex;P=C.enumeration.problem
    ids=[[sum(length,C.cells[1:k-1];init=0)+i for i in eachindex(layer)] for (k,layer) in enumerate(C.cells)]
    links=[Dict{Any,MorseVertex}() for layer in C.cells for c in layer]
    faces=morse_cell_faces(M;progress)
    # The real scalar centre acts trivially on every cell; retain one matrix
    # from each +/- pair for these orbit computations only, never in SL stabs.
    actions=[[collect(values(Dict(min(AV.matrix_key(AV.encode_matrix(P,s)),
        AV.matrix_key(AV.encode_matrix(P,-s)))=>s for s in S))) for S in layer] for layer in C.stabilizers]
    bases=Float64[]
    edge_orbits=NTuple{4,Int}[]
    negative_actions=[[[det(g)==-1 for g in S] for S in layer] for layer in actions]
    last_report=time()
    for (k,layer) in enumerate(C.cells),(i,c) in enumerate(layer)
        top=M.vertices[k][i]
        lower=links[ids[k][i]]
        for v in faces[k][i];lower[v.key]=v;end
        seen=Set{Any}()
        for v in faces[k][i]
            v.key in seen && continue
            # Equivalent incidences under Stab(top) give the same upper link
            # after expansion by Stab(v). Normalize only one from each orbit.
            flips=isodd(P.rank) # -I fixes every flag and has determinant -1 in odd rank.
            for (s,negative) in zip(actions[k][i],negative_actions[k][i])
                key=morse_act(M,s,v).key;push!(seen,key)
                negative && key==v.key && (flips=true)
            end
            push!(edge_orbits,(ids[k][i],ids[v.dim+1][v.orbit],v.shift,flips ? 1 : 2))
            upper=links[ids[v.dim+1][v.orbit]]
            inverse=inv(v.translate)
            for s in actions[v.dim+1][v.orbit]
                w=morse_act(M,s*inverse,top)
                upper[w.key]=w
            end
        end
        clear_action_cache!(M)
        F=sum(AV.act_form(P,c.point,inv(s)) for s in C.stabilizers[k][i])/length(C.stabilizers[k][i])
        push!(bases,morse_base_height(P,M.units,F))
        if progress && (time()-last_report>=10 || ids[k][i]==length(links))
            println("Morse link incidences: ",ids[k][i],"/",length(links));flush(stdout)
            last_report=time()
        end
    end
    (Dict("D"=>P.discriminant,"rank"=>P.rank,"link_model"=>"flag_graph",
        "base_heights"=>bases,"orbit_dimensions"=>[k-1 for (k,l) in enumerate(C.cells) for c in l],
        "flag_edge_orbits"=>edge_orbits,"links"=>MorseLinkGraphs(links,ids,M)),links,faces)
end

function flag_adjacency(vertices,edges)
    graph=Dict(v=>Set([v]) for v in vertices)
    for (u,v) in edges
        u==v && throw(ArgumentError("a link graph has a loop"))
        haskey(graph,u) && haskey(graph,v) || continue
        push!(graph[u],v);push!(graph[v],u)
    end
    graph
end

function collapse_flag_graph(vertices,edges;progress=false)
    graph=flag_adjacency(vertices,edges)
    steps=Union{NTuple{2,Int},NTuple{3,Int}}[]
    pending=sort!(collect(keys(graph)))
    edge_pending=BitSet(pending)
    last_report=time()
    while length(graph)>1
        # Removing a vertex can create new dominated vertices only among its
        # neighbours; new dominated edges have both endpoints there. Removing
        # an edge can create new dominated vertices or edges only at its two
        # endpoints. Revisit those neighbourhoods, not the entire graph.
        while !isempty(pending)
            next=BitSet()
            for v in pending
                haskey(graph,v) || continue
                for w in sort!(collect(graph[v]))
                    w==v && continue
                    issubset(graph[v],graph[w]) || continue
                    push!(steps,(v,w))
                    union!(next,graph[v]);union!(edge_pending,graph[v])
                    for u in graph[v];u==v || delete!(graph[u],v);end
                    delete!(graph,v)
                    break
                end
                if progress && time()-last_report>=10
                    println("Link collapse: ",length(graph)," vertices remain; ",length(steps)," steps");flush(stdout)
                    last_report=time()
                end
            end
            pending=collect(next)
        end
        length(graph)==1 && break
        isempty(edge_pending) && break
        u=first(edge_pending);delete!(edge_pending,u)
        haskey(graph,u) || continue
        next=BitSet()
        # Visit all neighbours: u may have been queued after an edge change
        # at a larger-numbered vertex. Process a whole neighbourhood before
        # revisiting dominated vertices, avoiding one global scan per edge.
        for v in sort!(collect(graph[u]))
            v==u && continue
            common=length(graph[u])<length(graph[v]) ? intersect(graph[u],graph[v]) : intersect(graph[v],graph[u])
            delete!(common,u);delete!(common,v)
            candidates=sort!(collect(common))
            w=findfirst(w->issubset(common,graph[w]),candidates)
            isnothing(w) && continue
            push!(steps,(u,v,candidates[w]))
            delete!(graph[u],v);delete!(graph[v],u)
            push!(next,u,v);push!(edge_pending,u,v)
        end
        pending=collect(next)
        if progress && time()-last_report>=10
            println("Link collapse: ",length(graph)," vertices remain; ",length(steps)," steps");flush(stdout)
            last_report=time()
        end
    end
    (;passed=length(graph)==1,steps,model=:graph,expanded_steps=length(steps),residual=sort!(collect(keys(graph))))
end

# Replay checks every cone condition afresh, independently of the search.
# For sigma={v} or {u,v}, pair every remaining simplex containing sigma but
# not the apex with its union with the apex, in descending dimension. This
# expands each recorded step into elementary simplicial collapses.
function replay_flag_collapse(vertices,edges,steps)
    graph=flag_adjacency(vertices,edges)
    for step in steps
        @assert length(step) in (2,3) && length(unique(step))==length(step)
        @assert all(v->haskey(graph,v),step)
        if length(step)==2
            v,w=step
            @assert w in graph[v]
            @assert all(u->u in graph[w],graph[v]) "vertex link is not the recorded cone"
            for u in collect(graph[v]);u==v || delete!(graph[u],v);end
            delete!(graph,v)
        else
            u,v,w=step
            @assert v in graph[u] && w in graph[u] && w in graph[v]
            @assert all(x->x==u || x==v || x in graph[w],intersect(graph[u],graph[v])) "edge link is not the recorded cone"
            delete!(graph[u],v);delete!(graph[v],u)
        end
    end
    @assert length(graph)==1 "collapse witness does not end at a point"
    true
end

# Translating every height by the same rational number preserves all link
# comparisons. Choose a regular level with the fewest SL vertex orbits. Each
# GL edge-flag orbit contributes its determinant-sign multiplicity whenever
# an integer translate crosses the level; changes occur only at the endpoint
# heights modulo one. This is an exact sweep, with no field-specific choices.
function least_vertex_level(bases,edges)
    events=zeros(BigInt,length(bases));count=big(0)
    for (top,lower,shift,weight) in edges
        a=bases[top];b=bases[lower]+shift
        low,high=a<b ? (top,lower) : (lower,top)
        lo,hi=minmax(a,b)
        count+=weight*(floor(BigInt,-lo)-floor(BigInt,-hi))
        events[low]+=weight;events[high]-=weight
    end
    before=count;best=count;offset=big(0)//big(1)
    fractions=mod.(bases,1);order=sortperm(fractions)
    for (j,i) in enumerate(order)
        count+=events[i]
        if count<best
            next=j==length(order) ? one(fractions[i]) : fractions[order[j+1]]
            offset=(fractions[i]+next)/2;best=count
        end
    end
    @assert count==before
    offset,before,best
end

function rational_morse_heights(data;bases=nothing)
    if isnothing(bases)
        scale=max(big(1_000_000),big(length(data["base_heights"]))+1)
        bases=[(round(BigInt,x*1000)*scale+i)//(1000scale) for (i,x) in enumerate(data["base_heights"])]
        if haskey(data,"flag_edge_orbits")
            offset,before,after=least_vertex_level(bases,data["flag_edge_orbits"])
            bases=bases .- offset
            data["level_selection"]=Dict("offset"=>string(offset),
                "zero_cells_before"=>before,"zero_cells_after"=>after)
        end
    end
    @assert length(data["links"])==length(bases) "incomplete collection of vertex links"
    @assert length(unique(mod.(bases,1)))==length(bases)
    @assert all(!isinteger,bases) "zero must be a regular level"
    bases
end

function ascending_vertices(data,bases,orbit,sgn;link=data["links"][orbit])
    delta=[bases[v["orbit"]]+v["shift"]-bases[orbit] for v in link["vertices"]]
    @assert all(!iszero,delta)
    findall(x->sgn*x>0,delta)
end

function certify_morse_graphs(data;bases=nothing,progress=false,on_certificate=nothing)
    bases=rational_morse_heights(data;bases)
    reports=Any[]
    last_report=time()
    for i in eachindex(data["links"])
        if progress && time()-last_report>=10
            println("Building Morse link graph: ",i,"/",length(data["links"]));flush(stdout)
            last_report=time()
        end
        link=data["links"][i]
        graph_vertices=length(link["vertices"]);graph_edges=length(link["edges"])
        for sgn in (-1,1)
            vertices=ascending_vertices(data,bases,i,sgn;link)
            if progress && time()-last_report>=10
                println("Certifying Morse link ",i,", direction ",sgn,": ",length(vertices)," vertices");flush(stdout)
                last_report=time()
            end
            result=if haskey(link,"symmetry_generators")
                orbit_collapse_flag_graph(vertices,link["edges"],link["symmetry_generators"],link["colours"];progress)
            else
                collapse_flag_graph(vertices,link["edges"];progress)
            end
            if !result.passed && haskey(link,"symmetry_generators")
                # A different greedy order can succeed where another stops.
                # Both searches use fresh graph data and require exact replay.
                result=collapse_flag_graph(vertices,link["edges"];progress)
            end
            @assert result.passed "No graph-collapse certificate for link $i, direction $sgn"
            @assert replay_morse_graph_certificate(link,vertices,result)
            report=(;orbit=i,sgn,result,graph_vertices,graph_edges)
            push!(reports,report)
            isnothing(on_certificate) || on_certificate(report,bases)
        end
        if progress && (i%100==0 || i==length(data["links"]) || time()-last_report>=10)
            println("Certified Morse links: ",i,"/",length(data["links"]));flush(stdout)
            last_report=time()
        end
    end
    (;passed=true,bases,reports)
end

# Dispatch on the recorded format, including after a plain-search fallback.
function replay_morse_graph_certificate(link,vertices,result)
    if result.model===:orbit
        replay_orbit_flag_collapse(vertices,link["edges"],link["symmetry_generators"],link["colours"],result.steps)
    elseif result.model===:graph
        replay_flag_collapse(vertices,link["edges"],result.steps)
    else
        throw(ArgumentError("unknown graph-collapse certificate model $(result.model)"))
    end
end
