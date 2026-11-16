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

function morse_link_graphs(M;progress=false)
    C=M.complex;P=C.enumeration.problem
    ids=[[sum(length,C.cells[1:k-1];init=0)+i for i in eachindex(layer)] for (k,layer) in enumerate(C.cells)]
    links=[Dict{Any,MorseVertex}() for layer in C.cells for c in layer]
    faces=morse_cell_faces(M;progress)
    for (k,layer) in enumerate(C.cells),(i,c) in enumerate(layer)
        top=M.vertices[k][i]
        lower=links[ids[k][i]]
        for v in faces[k][i]
            lower[v.key]=v
            upper=links[ids[v.dim+1][v.orbit]]
            for s in C.stabilizers[v.dim+1][v.orbit]
                w=morse_act(M,s*inv(v.translate),top)
                upper[w.key]=w
            end
        end
        clear_action_cache!(M)
    end
    data=Any[]
    bases=Float64[]
    for (k,layer) in enumerate(C.cells),(i,c) in enumerate(layer)
        vertices=sort!(collect(values(links[ids[k][i]]));by=v->v.key)
        # Incidence is reverse inclusion of complete minimum sets. An inverted
        # ray index finds supersets directly; testing every pair of vertices
        # would be quadratic in the size of a high-dimensional link.
        containing=Dict{Any,BitSet}()
        for (j,v) in enumerate(vertices),ray in v.key[2]
            push!(get!(containing,ray,BitSet()),j)
        end
        edges=Tuple{Int,Int}[]
        for (b,v) in enumerate(vertices)
            lists=sort!([containing[ray] for ray in v.key[2]];by=length)
            candidates=copy(first(lists))
            for list in lists[2:end];intersect!(candidates,list);end
            for a in candidates
                vertices[a].dim<v.dim && push!(edges,(a,b))
            end
        end
        push!(data,Dict("vertices"=>[Dict("orbit"=>ids[v.dim+1][v.orbit],"shift"=>v.shift) for v in vertices],"edges"=>edges))
        F=sum(AV.act_form(P,c.point,inv(s)) for s in C.stabilizers[k][i])/length(C.stabilizers[k][i])
        push!(bases,morse_base_height(P,M.units,F))
        if progress && (ids[k][i]%100==0 || ids[k][i]==length(links))
            println("Morse link graphs: ",ids[k][i],"/",length(links));flush(stdout)
        end
    end
    (Dict("D"=>P.discriminant,"rank"=>P.rank,"link_model"=>"flag_graph",
        "base_heights"=>bases,"orbit_dimensions"=>[k-1 for (k,l) in enumerate(C.cells) for c in l],
        "links"=>data),links,faces)
end

function flag_adjacency(vertices,edges)
    graph=Dict(v=>BitSet([v]) for v in vertices)
    for (u,v) in edges
        u==v && throw(ArgumentError("a link graph has a loop"))
        haskey(graph,u) && haskey(graph,v) || continue
        push!(graph[u],v);push!(graph[v],u)
    end
    graph
end

function collapse_flag_graph(vertices,edges)
    graph=flag_adjacency(vertices,edges)
    steps=Vector{Int}[]
    while length(graph)>1
        changed=false
        for v in sort!(collect(keys(graph)))
            haskey(graph,v) || continue
            for w in graph[v]
                w==v && continue
                issubset(graph[v],graph[w]) || continue
                push!(steps,[v,w])
                for u in graph[v];u==v || delete!(graph[u],v);end
                delete!(graph,v)
                changed=true
                break
            end
        end
        changed && continue
        for u in sort!(collect(keys(graph)))
            for v in graph[u]
                v>u || continue
                common=intersect(graph[u],graph[v])
                delete!(common,u);delete!(common,v)
                w=findfirst(w->issubset(common,graph[w]),collect(common))
                isnothing(w) && continue
                apex=collect(common)[w]
                push!(steps,[u,v,apex])
                delete!(graph[u],v);delete!(graph[v],u)
                changed=true
                break
            end
            changed && break
        end
        changed || break
    end
    (;passed=length(graph)==1,steps,residual=sort!(collect(keys(graph))))
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

function rational_morse_heights(data;bases=nothing)
    if isnothing(bases)
        scale=max(big(1_000_000),big(length(data["base_heights"]))+1)
        bases=[(round(BigInt,x*1000)*scale+i)//(1000scale) for (i,x) in enumerate(data["base_heights"])]
    end
    @assert length(data["links"])==length(bases) "incomplete collection of vertex links"
    @assert length(unique(mod.(bases,1)))==length(bases)
    @assert all(!isinteger,bases) "zero must be a regular level"
    bases
end

function ascending_vertices(data,bases,orbit,sgn)
    delta=[bases[v["orbit"]]+v["shift"]-bases[orbit] for v in data["links"][orbit]["vertices"]]
    @assert all(!iszero,delta)
    findall(x->sgn*x>0,delta)
end

function certify_morse_graphs(data;bases=nothing,progress=false)
    bases=rational_morse_heights(data;bases)
    reports=Any[]
    for (i,link) in enumerate(data["links"]),sgn in (-1,1)
        vertices=ascending_vertices(data,bases,i,sgn)
        result=collapse_flag_graph(vertices,link["edges"])
        @assert result.passed "No graph-collapse certificate for link $i, direction $sgn"
        @assert replay_flag_collapse(vertices,link["edges"],result.steps)
        push!(reports,(;orbit=i,sgn,result))
        if progress && (i%100==0 || i==length(data["links"])) && sgn==1
            println("Certified Morse links: ",i,"/",length(data["links"]));flush(stdout)
        end
    end
    (;passed=true,bases,reports)
end
