# Algebraic cancellation of isomorphic induced orientation-module summands.
# The unreduced object remains the geometric full or bounded-degree complex.
struct ReducedLevelComplex
    experiment::MorseExperiment
    certificate::Any
    flag_orbits::Vector{FlagOrbit}
    cells::Any
    stabilizers::Any
    boundaries::Any
    generators::Vector{Any}
    cancellation_pairs::Any
    original::Any
    original_indices::Vector{Vector{Int}}
end
level_orientation(L::ReducedLevelComplex,k,i,g)=level_orientation(L.original,k,L.original_indices[k][i],g)
level_model(L::ReducedLevelComplex)=level_model(L.original)
retained_parent_dimension(L::ReducedLevelComplex)=retained_parent_dimension(L.original)
function reduce_level_complex(L;progress=false)
    P=L.experiment.complex.enumeration.problem
    active=[trues(length(layer)) for layer in L.cells]
    boundaries=[[copy(b) for b in layer] for layer in L.boundaries]
    generators=unique([g for layer in L.stabilizers for S in layer for g in S] ∪
                      [g for layer in L.boundaries for b in layer for (j,g) in b])
    canonical_cache=Dict{Any,Any}()
    function canonical(k,i,g)
        key=(k,i,AV.matrix_key(AV.encode_matrix(P,g)))
        get!(canonical_cache,key) do
            matrices=[g*s for s in L.stabilizers[k][i]]
            keys=[AV.matrix_key(AV.encode_matrix(P,h)) for h in matrices]
            j=argmin(keys)
            (keys[j],matrices[j],level_orientation(L,k,i,L.stabilizers[k][i][j]))
        end
    end
    function simplify(k,terms)
        counts=Dict{Any,Int}();representatives=Dict{Any,Any}()
        for (j,g) in terms
            key,h,orientation=canonical(k,abs(j),g)
            id=(abs(j),key)
            counts[id]=get(counts,id,0)+sign(j)*orientation
            representatives[id]=h
        end
        result=Tuple{Int,Any}[]
        for (id,n) in counts
            for _ in 1:abs(n);push!(result,(sign(n)*id[1],representatives[id]));end
        end
        result
    end
    # Index all cells incident to each lower cell. Only changed rows and rows
    # whose pivot cost changes need reconsideration after a cancellation.
    incoming=[zeros(Int,length(layer)) for layer in L.cells]
    incident=[[Set{Int}() for _ in layer] for layer in L.cells]
    orders=[[length(S) for S in layer] for layer in L.stabilizers]
    for k in 2:length(L.cells),(i,b) in enumerate(boundaries[k]),(j,g) in b
        incoming[k-1][abs(j)]+=1
        push!(incident[k-1][abs(j)],i)
    end
    versions=[zeros(Int,length(layer)) for layer in L.cells]
    candidates=[Any[nothing for _ in layer] for layer in L.cells]
    heap=NTuple{4,Int}[]
    function enqueue!(entry)
        push!(heap,entry);i=length(heap)
        while i>1
            parent=i÷2
            heap[parent]<=entry && break
            heap[i]=heap[parent];i=parent
        end
        heap[i]=entry
    end
    function dequeue!()
        first=heap[1];last=pop!(heap)
        isempty(heap) && return first
        i=1
        while 2i<=length(heap)
            child=2i
            child<length(heap) && heap[child+1]<heap[child] && (child+=1)
            last<=heap[child] && break
            heap[i]=heap[child];i=child
        end
        heap[i]=last
        first
    end
    function refresh!(k,i)
        versions[k][i]+=1
        candidates[k][i]=nothing
        active[k][i] || return
        b=boundaries[k][i];counts=Dict{Int,Int}()
        for (j,g) in b;counts[abs(j)]=get(counts,abs(j),0)+1;end
        bestcost=typemax(Int)
        for (j,g) in b
            counts[abs(j)]==1 || continue
            orders[k][i]==orders[k-1][abs(j)] || continue
            @assert active[k-1][abs(j)] && incoming[k-1][abs(j)]>0
            cost=(incoming[k-1][abs(j)]-1)*(length(b)-1)
            if cost<bestcost
                candidates[k][i]=(j,g);bestcost=cost
            end
        end
        isnothing(candidates[k][i]) || enqueue!((bestcost,k,i,versions[k][i]))
    end
    function replace_boundary!(k,i,terms,affected)
        old=boundaries[k][i]
        columns=Set(abs(j) for (j,g) in old)
        union!(columns,(abs(j) for (j,g) in terms))
        for b in columns,q in incident[k-1][b];push!(affected,(k,q));end
        for (j,g) in old;incoming[k-1][abs(j)]-=1;end
        for b in columns;delete!(incident[k-1][b],i);end
        boundaries[k][i]=terms
        for (j,g) in terms
            incoming[k-1][abs(j)]+=1;push!(incident[k-1][abs(j)],i)
        end
        for b in columns,q in incident[k-1][b];push!(affected,(k,q));end
        push!(affected,(k,i))
    end
    for k in 2:length(L.cells),i in eachindex(L.cells[k]);refresh!(k,i);end
    cancellations=Any[]
    while !isempty(heap)
        cost,k,i,version=dequeue!()
        version==versions[k][i] && active[k][i] || continue
        j,g=candidates[k][i]
        b=abs(j);s=sign(j);inverse=inv(g)
        @assert Set(AV.matrix_key(AV.encode_matrix(P,g*h*inverse)) for h in L.stabilizers[k-1][b])==
                Set(AV.matrix_key(AV.encode_matrix(P,h)) for h in L.stabilizers[k][i])
        @assert all(level_orientation(L,k,i,h)==level_orientation(L,k-1,b,inverse*h*g)
                    for h in L.stabilizers[k][i]) "pivot orientation characters are incompatible"
        remainder=[(l,h) for (l,h) in boundaries[k][i] if abs(l)!=b]
        affected=Set{Tuple{Int,Int}}()
        for q in collect(incident[k-1][b])
            q==i && continue
            @assert active[k][q]
            hits=[(l,h) for (l,h) in boundaries[k][q] if abs(l)==b]
            terms=[(l,h) for (l,h) in boundaries[k][q] if abs(l)!=b]
            for (l,h) in hits,(m,t) in remainder
                push!(terms,(-sign(l)*s*m,h*inverse*t))
            end
            replace_boundary!(k,q,simplify(k-1,terms),affected)
        end
        if k<length(L.cells)
            for q in collect(incident[k][i])
                @assert active[k+1][q]
                terms=[term for term in boundaries[k+1][q] if abs(term[1])!=i]
                replace_boundary!(k+1,q,terms,affected)
            end
        end
        active[k][i]=false;active[k-1][b]=false
        replace_boundary!(k,i,Tuple{Int,Any}[],affected)
        # The cancelled lower row may itself have outgoing incidences.
        if k>2;replace_boundary!(k-1,b,Tuple{Int,Any}[],affected);end
        for (degree,row) in affected;refresh!(degree,row);end
        push!(cancellations,(k-1,i,b,g))
        if progress && length(cancellations)%100==0
            println("Equivariant chain cancellations: ",length(cancellations),"; remaining ",count.(active));flush(stdout)
        end
    end
    remaps=[Dict(old=>new for (new,old) in enumerate(findall(a))) for a in active]
    cells=[layer[active[k]] for (k,layer) in enumerate(L.cells)]
    stabs=[layer[active[k]] for (k,layer) in enumerate(L.stabilizers)]
    result=[[[(sign(j)*remaps[k-1][abs(j)],g) for (j,g) in b] for (i,b) in enumerate(layer) if active[k][i]]
        for (k,layer) in enumerate(boundaries)]
    if progress;println("Reduced level complex: ",length.(cells));flush(stdout);end
    ReducedLevelComplex(L.experiment,L.certificate,L.flag_orbits,cells,stabs,result,generators,cancellations,L,findall.(active))
end
