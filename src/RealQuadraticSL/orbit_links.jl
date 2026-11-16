# A cell stabilizer preserves the dimension colours of its link. A dominated
# vertex and its apex have different colours; a dominated edge and its apex
# have three different colours. Consequently deleting one orbit member cannot
# delete an apex or a cone edge needed by another member of the same orbit.
# Expand each orbit removal into ordinary steps, checking every cone condition;
# the separate replay in flag_links.jl checks the final witness independently.
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
    steps=Union{NTuple{2,Int},NTuple{3,Int}}[]
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
                for w in sort!(collect(graph[v]);by=w->(-length(graph[w]),w))
                    w==v && continue
                    issubset(graph[v],graph[w]) || continue
                    @assert colours[v]!=colours[w]
                    pairs=images((v,w));@assert length(pairs)==length(members[o])
                    affected=BitSet()
                    for (x,y) in pairs
                        @assert haskey(graph,x) && haskey(graph,y) && issubset(graph[x],graph[y])
                        push!(steps,(x,y))
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
                    println("Orbit collapse: ",length(graph)," vertices; ",length(steps)," steps");flush(stdout);last_report=time()
                end
            end
            pending=collect(next)
        end
        length(graph)==1 && break
        if !initialized;foreach(enqueue,eachindex(reps));initialized=true;end
        isempty(degrees) && break
        while !haskey(buckets,minimum_degree[]) || isempty(buckets[minimum_degree[]]);minimum_degree[]+=1;end
        o=first(buckets[minimum_degree[]]);dequeue(o);u=reps[o]
        next=BitSet()
        for v in sort!(collect(graph[u]);by=v->(length(graph[v]),v))
            (v==u || !(v in graph[u])) && continue
            common=length(graph[u])<length(graph[v]) ? intersect(graph[u],graph[v]) : intersect(graph[v],graph[u])
            delete!(common,u);delete!(common,v)
            candidates=sort!(collect(common);by=w->(-length(graph[w]),w))
            w=findfirst(w->issubset(common,graph[w]),candidates)
            isnothing(w) && continue
            apex=candidates[w]
            @assert length(unique((colours[u],colours[v],colours[apex])))==3
            for (x,y,z) in images((u,v,apex))
                @assert y in graph[x] && z in graph[x] && z in graph[y]
                @assert all(t->t==x || t==y || t in graph[z],intersect(graph[x],graph[y]))
                push!(steps,(x,y,z));delete!(graph[x],y);delete!(graph[y],x)
            end
            push!(next,o,orbit[v]);enqueue(o);enqueue(orbit[v])
        end
        pending=collect(next)
        if progress && time()-last_report>=10
            println("Orbit collapse: ",length(graph)," vertices; ",length(steps)," steps");flush(stdout);last_report=time()
        end
    end
    (;passed=length(graph)==1,steps,residual=sort!(collect(keys(graph))))
end
