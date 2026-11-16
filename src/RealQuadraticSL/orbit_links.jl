# In a properly coloured graph, a cone apex cannot share its colour with
# another vertex of the link. This necessary condition leaves at most one
# candidate per cell dimension, even when a common neighbourhood is huge.
function unique_colour_candidates(vertices,colours)
    singletons=Dict{Int,Int}()
    for v in vertices
        c=colours[v]
        singletons[c]=haskey(singletons,c) ? 0 : v
    end
    [v for v in values(singletons) if v!=0]
end

# Proper dimension colours are invariant under the stabilizer. Two distinct
# members of a vertex or edge orbit cannot occur in a common simplex. Their
# removals therefore commute: deleting one cannot change the link of another.
# Record one witness per orbit and independently verify the graph action and
# every complete removal orbit. A cone link gives an immediate witness. If a
# vertex link instead collapses to a point, coning its collapse pairs with the
# removed vertex collapses its star relative to its link, ending with one edge.
function orbit_collapse_flag_graph(vertices,edges,generators,colours;progress=false)
    graph=flag_adjacency(vertices,edges)
    @assert all(colours[u]!=colours[v] for (u,v) in edges)
    @assert all(sort(g)==collect(eachindex(colours)) for g in generators)
    @assert all(colours[g[v]]==colours[v] for g in generators for v in eachindex(colours))
    orbit=zeros(Int,length(colours));members=Vector{Int}[]
    for v in sort!(collect(keys(graph)))
        orbit[v]>0 && continue
        push!(members,[v]);o=length(members);orbit[v]=o
        for x in members[o],g in generators
            y=g[x];@assert haskey(graph,y)
            if orbit[y]==0;orbit[y]=o;push!(members[o],y);end
            @assert orbit[y]==o
        end
    end
    reps=first.(members)
    if progress;println("Symmetry orbits: ",length(members)," / ",length(graph));flush(stdout);end
    function images(tuple)
        seen=Set([tuple[1:end-1]]);queue=[tuple]
        for item in queue,g in generators
            image=map(v->g[v],item)
            key=image[1:end-1]
            key in seen && continue
            push!(seen,key);push!(queue,image)
        end
        queue
    end
    orbit_steps=Any[]
    expanded_steps=0
    pending=collect(eachindex(reps))
    degrees=Dict{Int,Int}();buckets=Dict{Int,Set{Int}}();minimum_degree=Ref(typemax(Int))
    initialized=false
    function dequeue(o)
        old=pop!(degrees,o,0);old>0 && delete!(buckets[old],o)
    end
    function enqueue(o)
        v=reps[o]
        haskey(graph,v) || (dequeue(o);return)
        degree=length(graph[v]);old=get(degrees,o,0)
        old==degree && return
        old>0 && delete!(buckets[old],o)
        push!(get!(Set{Int},buckets,degree),o);degrees[o]=degree
        minimum_degree[]=min(minimum_degree[],degree)
    end
    last_report=time()
    while length(graph)>1
        while !isempty(pending)
            next=BitSet()
            for o in pending
                v=reps[o];haskey(graph,v) || continue
                for w in sort!(unique_colour_candidates(graph[v],colours);by=w->(-length(graph[w]),w))
                    w==v && continue
                    issubset(graph[v],graph[w]) || continue
                    @assert colours[v]!=colours[w]
                    pairs=images((v,w));@assert length(pairs)==length(members[o])
                    push!(orbit_steps,(v,w))
                    affected=BitSet()
                    for (x,y) in pairs
                        @assert haskey(graph,x) && haskey(graph,y) && y in graph[x]
                        expanded_steps+=1
                        for u in graph[x]
                            u==x && continue
                            push!(affected,orbit[u]);delete!(graph[u],x)
                        end
                        delete!(graph,x)
                    end
                    dequeue(o);union!(next,affected)
                    initialized && foreach(enqueue,affected)
                    break
                end
                if progress && time()-last_report>=10
                    println("Orbit collapse: ",length(graph)," vertices; ",expanded_steps," steps");flush(stdout);last_report=time()
                end
            end
            pending=collect(next)
        end
        length(graph)==1 && break
        if !initialized;foreach(enqueue,eachindex(reps));initialized=true;end
        if isempty(degrees)
            # A collapsible link suffices even when it has no cone vertex.
            candidates=sort!([o for o in eachindex(reps) if haskey(graph,reps[o])];by=o->(length(graph[reps[o]]),o))
            removed=false
            for o in candidates
                v=reps[o];neighbours=setdiff(graph[v],Set([v]))
                link_edges=[(x,y) for x in neighbours for y in intersect(neighbours,graph[x]) if x<y]
                trial=collapse_flag_graph(neighbours,link_edges)
                trial.passed || continue
                push!(orbit_steps,Dict("vertex"=>v,"link_steps"=>trial.steps))
                affected=BitSet()
                for x in members[o]
                    for u in graph[x]
                        u==x && continue
                        push!(affected,orbit[u]);delete!(graph[u],x)
                    end
                    delete!(graph,x);expanded_steps+=1
                end
                foreach(enqueue,affected);pending=collect(affected);removed=true
                if progress;println("Collapsible-link orbit: ",length(members[o])," vertices removed; ",length(graph)," remain; ",length(trial.steps)," inner steps");flush(stdout);end
                break
            end
            removed || break
            continue
        end
        while !haskey(buckets,minimum_degree[]) || isempty(buckets[minimum_degree[]]);minimum_degree[]+=1;end
        o=first(buckets[minimum_degree[]]);dequeue(o);u=reps[o]
        next=BitSet()
        for v in sort!(collect(graph[u]);by=v->(length(graph[v]),v))
            (v==u || !(v in graph[u])) && continue
            common=length(graph[u])<length(graph[v]) ? intersect(graph[u],graph[v]) : intersect(graph[v],graph[u])
            delete!(common,u);delete!(common,v)
            candidates=sort!(unique_colour_candidates(common,colours);by=w->(-length(graph[w]),w))
            w=findfirst(w->issubset(common,graph[w]),candidates)
            isnothing(w) && continue
            apex=candidates[w]
            @assert length(unique((colours[u],colours[v],colours[apex])))==3
            push!(orbit_steps,(u,v,apex))
            for (x,y,z) in images((u,v,apex))
                @assert y in graph[x] && z in graph[x] && z in graph[y]
                expanded_steps+=1;delete!(graph[x],y);delete!(graph[y],x)
            end
            push!(next,o,orbit[v]);enqueue(o);enqueue(orbit[v])
        end
        pending=collect(next)
        if progress && time()-last_report>=10
            println("Orbit collapse: ",length(graph)," vertices; ",expanded_steps," steps");flush(stdout);last_report=time()
        end
    end
    (;passed=length(graph)==1,steps=orbit_steps,expanded_steps,model=:orbit,residual=sort!(collect(keys(graph))))
end

# Replay validates the proposed group action against the graph itself. It
# rebuilds each removal orbit by a separate depth-first traversal, checks the
# representative cone, and uses invariant proper colours to ensure that orbit
# members cannot delete one another's apices or cone edges.
function replay_orbit_flag_reduction(vertices,edges,generators,colours,steps)
    G=flag_adjacency(vertices,edges)
    @assert all(colours[u]!=colours[v] for (u,v) in edges)
    for permutation in generators
        @assert sort(permutation)==collect(eachindex(colours))
        @assert all(colours[permutation[v]]==colours[v] for v in eachindex(colours))
        @assert all(haskey(G,permutation[v]) for v in keys(G))
        for (u,neighbours) in G
            image=G[permutation[u]]
            @assert all(permutation[v] in image for v in neighbours) "not a graph automorphism"
        end
    end
    for entry in steps
        if entry isa AbstractDict
            @assert Set(keys(entry))==Set(["vertex","link_steps"])
            v=Int(entry["vertex"]);@assert haskey(G,v)
            neighbours=setdiff(G[v],Set([v]))
            link_edges=[(x,y) for x in neighbours for y in intersect(neighbours,G[x]) if x<y]
            @assert replay_flag_collapse(neighbours,link_edges,entry["link_steps"])
            stack=[v];seen=Set{Int}()
            while !isempty(stack)
                x=pop!(stack);x in seen && continue;push!(seen,x)
                for g in generators;push!(stack,g[x]);end
            end
            for x in seen
                @assert haskey(G,x)
                for u in G[x];u==x || delete!(G[u],x);end
                delete!(G,x)
            end
            continue
        end
        step=Tuple(Int.(entry))
        @assert length(step) in (2,3)
        @assert length(unique(colours[v] for v in step))==length(step)
        @assert all(v->haskey(G,v),step)
        if length(step)==2
            v,w=step
            @assert w in G[v] && all(x->x in G[w],G[v])
        else
            u,v,w=step
            @assert v in G[u] && w in G[u] && w in G[v]
            @assert all(x->x==u || x==v || x in G[w],intersect(G[u],G[v]))
        end
        stack=[step];seen=Set{Tuple}();images=Tuple[]
        while !isempty(stack)
            item=pop!(stack);face=item[1:end-1]
            face in seen && continue
            push!(seen,face);push!(images,item)
            for permutation in generators
                push!(stack,map(v->permutation[v],item))
            end
        end
        if length(step)==2
            for (v,w) in images
                @assert haskey(G,v) && haskey(G,w) && w in G[v]
                for x in G[v];x==v || delete!(G[x],v);end
                delete!(G,v)
            end
        else
            for (u,v,w) in images
                @assert v in G[u] && w in G[u] && w in G[v]
                delete!(G[u],v);delete!(G[v],u)
            end
        end
    end
    G
end

function replay_orbit_flag_collapse(vertices,edges,generators,colours,steps)
    residual=replay_orbit_flag_reduction(vertices,edges,generators,colours,steps)
    @assert length(residual)==1
    true
end

# The residual graph need not be a comparability graph after edge collapses.
# Compute its clique dimension directly, pruning by the invariant proper colours.
# This is an exact bound on the remaining simplicial core, not a homology test.
function flag_graph_dimension(graph,colours)
    @assert all(u==v || colours[u]!=colours[v] for (u,neighbours) in graph for v in neighbours)
    best=Ref(0)
    function search(candidates,depth)
        best[]=max(best[],depth)
        depth+length(candidates)<=best[] && return
        depth+length(unique(colours[v] for v in candidates))<=best[] && return
        while !isempty(candidates)
            v=pop!(candidates)
            search([w for w in candidates if w in graph[v]],depth+1)
            depth+length(candidates)<=best[] && return
        end
    end
    search(sort!(collect(keys(graph));by=v->length(graph[v])),0)
    best[]-1
end
