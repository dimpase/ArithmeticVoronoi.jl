# Exact traversal of the relative Ryshkov polyhedron, using Julia/Hecke/Oscar
# objects throughout. See Schuermann, "Enumerating perfect forms", section 7
# (https://arxiv.org/abs/0901.1587) for the T-perfect-form framework.

function native_direction(P,coefficients)
    sum(c*B for (c,B) in zip(coefficients,P.form_basis))
end

native_has_negative_direction(A)=last(Hecke.signature_tuple(quadratic_space(QQ,A;cached=false)))>0

"""Move to the first new minimum along an exact rational form-space direction."""
function native_flip(A,R,data=shortest_data(A))
    values=[quadratic_value(R,v) for v in data.vectors]
    ensure(all(>=(0),values) && any(iszero,values),"flip must preserve a nonempty minimum face")
    ensure(native_has_negative_direction(R),"flip direction has no negative eigenvalue")
    m=data.minimum
    lower,upper=QQ(0),QQ(1)
    local H,next
    # Find a positive-definite point beyond the first new minimum. Keeping a
    # lower endpoint inside the Ryshkov polyhedron avoids crossing its PD boundary.
    while true
        H=A+upper*R
        if !positive_definite(H)
            upper=(lower+upper)/2
            continue
        end
        next=shortest_data(H)
        next.minimum<m && break
        ensure(next.minimum==m,"flip lost its preserved minimum face")
        # A trial point can land exactly on the new vertex.
        any(v->quadratic_value(R,v)<0,next.vectors) && return (;gram=H,data=next)
        lower=upper
        upper*=2
    end
    # Each vector shorter than m gives an exact upper bound on the first hit.
    # These bounds decrease inside a compact PD segment, so only finitely many
    # lattice vectors can occur. No floating-point tolerance or deadline is used.
    while next.minimum<m
        bound=minimum((m-quadratic_value(A,v))/quadratic_value(R,v) for v in next.vectors)
        ensure(0<bound<upper,"new short vector did not shorten the flip interval")
        upper=bound
        H=A+upper*R
        next=shortest_data(H)
    end
    ensure(next.minimum==m,"flip changed the arithmetic minimum")
    (;gram=H,data=next)
end

function native_initial_form(P;progress=false)
    A=trace_gram(P,identity_matrix(P.field,P.rank))
    data=shortest_data(A)
    dimension=length(P.form_basis)
    while true
        rays=evaluation_rays(P,data.vectors)
        dimension-rank(rays)==0 && return (;gram=A/data.minimum,data=shortest_data(A/data.minimum))
        nullity,N=nullspace(rays)
        ensure(nullity>0,"imperfect form has no deformation direction")
        R=native_direction(P,[N[i,1] for i in 1:dimension])
        native_has_negative_direction(R) || (R=-R)
        next=native_flip(A,R,data)
        ensure(rank(evaluation_rays(P,next.data.vectors))>rank(rays),"initial-form deformation did not increase perfection rank")
        A,data=next.gram,next.data
        progress && println("Native perfect forms: seed rank ",dimension-nullity," -> ",rank(evaluation_rays(P,data.vectors)),"/",dimension)
        progress && flush(stdout)
    end
end

function native_vertex(P,A,data)
    ensure(data.minimum==1,"native representatives must have minimum one")
    F=form_from_trace(P,A)
    V=data.vectors
    rays=evaluation_rays(P,V)
    ensure(rank(rays)==length(P.form_basis),"native form is not perfect in its trace-form space")
    stabilizer=finite_stabilizer(P,F)
    closure=finite_matrix_closure([encode_matrix(P,g) for g in stabilizer.generators],P.rank*P.degree,stabilizer.order)
    ensure(haskey(closure.seen,matrix_key(-identity_matrix(QQ,P.rank*P.degree))),"native stabilizer lost central -I")
    index=Dict(Tuple(v)=>i for (i,v) in enumerate(V))
    permutations=[[index[signed_key(apply_integral(M,v))] for v in V] for M in closure.elements]
    all_facets=sort!(exact_facets(P,rays);by=f->Tuple(f.indices))
    facet_lookup=Dict(Tuple(f.indices)=>i for (i,f) in enumerate(all_facets))
    labels=Dict{Tuple,Tuple{Int,Int}}()
    orbit_facets=Int[]
    for (i,f) in enumerate(all_facets)
        haskey(labels,Tuple(f.indices)) && continue
        push!(orbit_facets,i)
        orbit=length(orbit_facets)
        for (h,permutation) in enumerate(permutations)
            image=Tuple(sort(permutation[f.indices]))
            ensure(haskey(facet_lookup,image),"native stabilizer does not permute exact facets")
            if haskey(labels,image)
                ensure(first(labels[image])==orbit,"native facet orbits overlap")
            else
                labels[image]=(orbit,h)
            end
        end
    end
    labelled=Any[]
    for f in all_facets
        orbit,h=labels[Tuple(f.indices)]
        permutation=invperm(permutations[h])
        ensure(sort(permutation[f.indices])==all_facets[orbit_facets[orbit]].indices,"native facet reduction has the wrong direction")
        push!(labelled,(;indices=f.indices,normal=f.normal,neighbour_orbit=orbit,
            to_orbit_representative=decode_matrix(P,inv(closure.elements[h])),
            stabilizer_word=-reverse(closure.words[h]),
            induced_minimal_vector_permutation=permutation))
    end
    record=PerfectFormRecord(F,A,QQ(1),V,rays,Any[stabilizer.generators...],stabilizer.order,labelled)
    (;record,orbit_facets,data)
end

# Determinant and kissing number are invariant under the allowed integral
# isometries. They only filter candidates; Hecke decides K-linear equivalence.
native_invariant(A,data)=(det(A),length(data.vectors))

function native_neighbour(P,source,target,facet,H,data,g,i,j)
    ensure(data.minimum==1 && positive_definite(H),"native neighbour has an invalid minimum or positivity")
    ensure(H!=source.trace_gram && rank(evaluation_rays(P,data.vectors))==length(P.form_basis),"native neighbour is not a distinct perfect form")
    shared=findall(v->quadratic_value(H,v)==1,source.minimal_vectors)
    ensure(shared==facet.indices,"native neighbour does not share exactly the crossed facet")
    HF=form_from_trace(P,H)
    direction=form_coordinates(P,HF-source.form)
    pivot=findfirst(!iszero,facet.normal)
    scale=direction[pivot]/facet.normal[pivot]
    ensure(scale>0 && direction==scale*facet.normal,"native neighbour is not on the outward facet ray")
    ensure(act_form(P,target.form,g)==HF,"native transporter has the wrong congruence direction")
    R=encode_matrix(P,g)
    check_group_matrix(P,R)
    index=Dict(Tuple(v)=>k for (k,v) in enumerate(target.minimal_vectors))
    ensure(Set(signed_key(apply_integral(R,v)) for v in data.vectors)==Set(keys(index)),"native transporter lost a minimal vector")
    mapping=Int[]
    for k in facet.indices
        image=apply_integral(R,source.minimal_vectors[k])
        key=signed_key(image)
        push!(mapping,Tuple(image)==key ? index[key] : -index[key])
    end
    ensure(sort(abs.(mapping)) in [f.indices for f in target.facets],"native transporter does not map the common facet to a target facet")
    NeighbourRecord(i,facet.indices,j,H,g,mapping)
end

function enumerate_native(P::ArithmeticVoronoiProblem;progress=false)
    start=time_ns()
    progress && println("Native perfect forms: constructing a starting form for discriminant ",P.discriminant)
    progress && flush(stdout)
    seed=native_initial_form(P;progress)
    vertices=[native_vertex(P,seed.gram,seed.data)]
    buckets=Dict(native_invariant(seed.gram,seed.data)=>[1])
    neighbours=NeighbourRecord[]
    edges_by_source=Vector{NeighbourRecord}[]
    coverage=Any[]
    cursor=1
    while cursor<=length(vertices)
        vertex=vertices[cursor]
        source=vertex.record
        edges=NeighbourRecord[]
        for facet_id in vertex.orbit_facets
            facet=source.facets[facet_id]
            next=native_flip(source.trace_gram,native_direction(P,facet.normal),vertex.data)
            H=form_from_trace(P,next.gram)
            key=native_invariant(next.gram,next.data)
            candidates=get!(buckets,key,Int[])
            target=0
            transporter=nothing
            for j in candidates
                g=equivalence(P,vertices[j].record.form,H)
                isnothing(g) && continue
                target,transporter=j,g
                break
            end
            if target==0
                push!(vertices,native_vertex(P,next.gram,next.data))
                target=length(vertices)
                push!(candidates,target)
                transporter=identity_matrix(P.field,P.rank)
            end
            edge=native_neighbour(P,source,vertices[target].record,facet,next.gram,next.data,transporter,cursor,target)
            push!(edges,edge)
            push!(neighbours,edge)
            if progress
                println("Native perfect forms: orbit ",cursor,", facet orbit ",length(edges),"/",length(vertex.orbit_facets),
                    "; ",length(vertices)," forms found, ",length(vertices)-cursor," awaiting traversal")
                flush(stdout)
            end
        end
        push!(edges_by_source,edges)
        push!(coverage,Dict("representative"=>cursor,"facets"=>length(source.facets),"facet_orbits"=>length(edges),"unresolved_facets"=>0))
        cursor+=1
    end
    # All facet orbits have been visited. Check each reverse incidence in the
    # target's own facet labelling before declaring the graph closed.
    for edge in neighbours
        target=vertices[edge.target].record
        facet=only(f for f in target.facets if f.indices==sort(abs.(edge.facet_map)))
        ensure(edges_by_source[edge.target][facet.neighbour_orbit].target==edge.source,"native reverse facet identification does not close")
    end
    provenance=Dict{String,Any}("backend"=>"native","algorithm"=>"relative Ryshkov neighbour traversal v1",
        "schema"=>SCHEMA,"from_cache"=>false,"backend_seconds"=>(time_ns()-start)/1e9,
        "julia_version"=>string(VERSION),"oscar_version"=>string(pkgversion(Oscar)),"hecke_version"=>string(pkgversion(Hecke)),
        "group_constraints"=>P.group_constraints,"action_convention"=>P.action_convention,
        "minimum_and_minimal_vectors"=>"exhaustive Hecke shortest vectors; exact squared trace lengths",
        "stabilizers"=>"Hecke K-linear automorphism groups; generated closure and central -I checked",
        "facets"=>"exact OSCAR/polymake cone facets; every facet covered by one stabilizer orbit",
        "neighbours"=>"exact minimum-preserving flips, common faces and full minimum-set transport; reverse incidences checked",
        "closure"=>coverage,"global_argument"=>"connectedness of the arithmetic perfect-form neighbour graph; all facet orbits exhausted under the specified K-linear GL group",
        "presentation_status"=>"not constructed; no ordered two-cell attaching loops yet")
    ArithmeticVoronoiEnumeration(P,[v.record for v in vertices],neighbours,Any[],:complete,provenance)
end
