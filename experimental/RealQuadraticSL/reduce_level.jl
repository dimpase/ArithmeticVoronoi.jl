# Algebraic cancellation of isomorphic permutation-module summands.
# The unreduced object remains the actual certified contractible cell complex.
struct ReducedLevelComplex
    experiment::MorseExperiment
    certificate::Any
    flag_orbits::Vector{FlagOrbit}
    cells::Any
    stabilizers::Any
    boundaries::Any
    generators::Vector{Any}
    cancellation_pairs::Any
end
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
            (keys[j],matrices[j])
        end
    end
    function simplify(k,terms)
        counts=Dict{Any,Int}();representatives=Dict{Any,Any}()
        for (j,g) in terms
            key,h=canonical(k,abs(j),g)
            id=(abs(j),key)
            counts[id]=get(counts,id,0)+sign(j)
            representatives[id]=h
        end
        result=Tuple{Int,Any}[]
        for (id,n) in counts
            for _ in 1:abs(n);push!(result,(sign(n)*id[1],representatives[id]));end
        end
        result
    end
    cancellations=Any[]
    while true
        best=nothing;bestcost=typemax(Int)
        for k in 2:length(L.cells)
            incoming=zeros(Int,length(L.cells[k-1]))
            for (i,b) in enumerate(boundaries[k])
                active[k][i] || continue
                for (j,g) in b;incoming[abs(j)]+=1;end
            end
            for (i,b) in enumerate(boundaries[k])
                active[k][i] || continue
                counts=Dict{Int,Int}()
                for (j,g) in b;counts[abs(j)]=get(counts,abs(j),0)+1;end
                for (j,g) in b
                    counts[abs(j)]==1 || continue
                    length(L.stabilizers[k][i])==length(L.stabilizers[k-1][abs(j)]) || continue
                    cost=(incoming[abs(j)]-1)*(length(b)-1)
                    if cost<bestcost
                        best=(k,i,j,g);bestcost=cost
                    end
                end
            end
        end
        isnothing(best) && break
        k,i,j,g=best
        b=abs(j);s=sign(j);inverse=inv(g)
        @assert Set(AV.matrix_key(AV.encode_matrix(P,inverse*h*g)) for h in L.stabilizers[k-1][b])==
                Set(AV.matrix_key(AV.encode_matrix(P,h)) for h in L.stabilizers[k][i])
        remainder=[(l,h) for (l,h) in boundaries[k][i] if abs(l)!=b]
        for q in eachindex(boundaries[k])
            q==i && continue
            active[k][q] || continue
            hits=[(l,h) for (l,h) in boundaries[k][q] if abs(l)==b]
            isempty(hits) && continue
            terms=[(l,h) for (l,h) in boundaries[k][q] if abs(l)!=b]
            for (l,h) in hits,(m,t) in remainder
                push!(terms,(-sign(l)*s*m,h*inverse*t))
            end
            boundaries[k][q]=simplify(k-1,terms)
        end
        if k<length(L.cells)
            for q in eachindex(boundaries[k+1])
                active[k+1][q] || continue
                filter!(x->abs(x[1])!=i,boundaries[k+1][q])
            end
        end
        active[k][i]=false;active[k-1][b]=false
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
    ReducedLevelComplex(L.experiment,L.certificate,L.flag_orbits,cells,stabs,result,generators,cancellations)
end
