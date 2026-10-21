export validate_export, ArithmeticValidationError

struct ArithmeticValidationError <: Exception
    reason::String
end
Base.showerror(io::IO,e::ArithmeticValidationError)=print(io,e.reason)
ensure(test,message)=test || throw(ArithmeticValidationError(message))

function rational_wire(s)
    s isa AbstractString && occursin(r"^-?[0-9]+(/[1-9][0-9]*)?$",s) || throw(ArithmeticValidationError("expected an exact integer/rational string"))
    parts=split(s,'/');a=ZZ(parts[1]);b=length(parts)==1 ? ZZ(1) : ZZ(parts[2])
    QQ(a,b)
end
function matrix_wire(rows)
    ensure(rows isa AbstractVector && !isempty(rows),"empty or malformed matrix")
    n=length(first(rows));ensure(n>0 && all(length(r)==n for r in rows),"ragged matrix")
    matrix(QQ,length(rows),n,[rational_wire(rows[i][j]) for i in eachindex(rows) for j in 1:n])
end
matrix_rows(A)=[[A[i,j] for j in 1:ncols(A)] for i in 1:nrows(A)]
matrix_key(A)=Tuple(BigInt(ZZ(x)) for x in A)
function form_coordinates(P,F)
    P.degree==1 && return QQFieldElem[F[1,1],F[1,2],F[2,2]]
    is_hermitian(P) && return QQFieldElem[coeff(F[1,1],0),coeff(F[1,2],0),coeff(F[1,2],1),coeff(F[2,2],0)]
    QQFieldElem[coeff(F[i,j],k) for (i,j) in [(1,1),(1,2),(2,2)] for k in 0:1]
end

function finite_matrix_closure(generators,dimension,expected_order)
    identity=identity_matrix(QQ,dimension)
    elements=[identity];words=[Int[]];seen=Dict(matrix_key(identity)=>1)
    cursor=1
    while cursor<=length(elements)
        for (j,g) in enumerate(generators)
            next=elements[cursor]*g;key=matrix_key(next)
            if !haskey(seen,key)
                ensure(length(elements)<expected_order,"generated stabilizer exceeds its claimed order")
                push!(elements,next);push!(words,[words[cursor];j]);seen[key]=length(elements)
            end
        end
        cursor+=1
    end
    ensure(length(elements)==expected_order,"stabilizer generators do not generate the claimed order")
    (;elements,words,seen)
end

function exact_facets(P,rays)
    cone=positive_hull(rays)
    ensure(dim(cone)==length(P.form_basis),"Voronoi cone is not full-dimensional in its form space")
    inequalities=matrix(QQ,Oscar.linear_inequality_matrix(facets(cone)))
    result=Any[]
    for j in 1:nrows(inequalities)
        normal=[inequalities[j,k] for k in 1:ncols(inequalities)]
        values=[sum(normal[k]*rays[i,k] for k in eachindex(normal)) for i in 1:nrows(rays)]
        if all(v<=0 for v in values);normal=-normal;values=-values;end
        ensure(all(v>=0 for v in values) && any(v>0 for v in values),"invalid exact facet normal")
        indices=findall(iszero,values)
        ensure(rank(rays[indices,:])==length(P.form_basis)-1,"facet has incorrect rank")
        push!(result,(;indices,normal))
    end
    result
end

function checked_vectors(P,raw,gram,minimum)
    data=shortest_data(gram)
    ensure(data.minimum==minimum,"backend minimum disagrees with exhaustive Hecke shortest vectors")
    M=matrix_wire(raw)
    ensure(all(isone(denominator(x)) for x in M),"nonintegral minimal vector")
    supplied=[internal_vector(P,v) for v in matrix_rows(M)]
    expected=Set(Tuple(v) for v in [data.vectors;-data.vectors])
    ensure(Set(Tuple(v) for v in supplied)==expected,"backend minimal vectors are not the exhaustive signed Hecke list")
    (;vectors=data.vectors,supplied)
end

"""Independent exact validation of the exported closure certificate, not a process-exit check."""
function validate_export(P::ArithmeticVoronoiProblem,data)
    ensure(get(data,"schema",nothing)=="av_polyhedral_export_v1","unsupported export schema")
    ensure(get(data,"completed",false)===true && isempty(data["pending_facets"]),"incomplete backend enumeration")
    ensure(data["arithmetic"]=="gmp","unexpected arithmetic")
    ensure(data["boundary_square"]=="not_computed","this importer accepts perfect-form data only; invalid or unsupported boundary assertion")
    ensure(data["action"]=="H=M*F*transpose(M); row vectors act by right multiplication","unknown backend action convention")
    ensure(data["minimal_vector_convention"]=="both signs in conversion_and_duplication order","unknown signed-vector convention")
    basis=[internal_gram(P,matrix_wire(B)) for B in data["form_basis"]]
    ensure(length(basis)==length(P.form_basis),"wrong form-space dimension")
    for B in basis;form_from_trace(P,B);end
    positions=[(i,j) for i in 1:P.rank*P.degree for j in 1:P.rank*P.degree]
    ensure(rank(matrix(QQ,[B[i,j] for B in basis,(i,j) in positions]))==length(basis),"backend form basis is dependent")
    if P.degree==2
        comms=data["commuting_matrices"]
        ensure(length(comms)==1 && internal_gram(P,transpose(matrix_wire(only(comms))))==P.multiplication,"wrong field-linearity constraint or multiplication convention")
    else
        ensure(isempty(data["commuting_matrices"]),"unexpected classical commuting constraints")
    end
    raw=data["representatives"];ensure(!isempty(raw),"no starting perfect form")
    records=PerfectFormRecord[];actions=Any[];signed_index_maps=Any[];facets_all=Any[]
    for (i,entry) in enumerate(raw)
        ensure(entry["id"]==i,"representative IDs must be consecutive and one-based")
        A=internal_gram(P,matrix_wire(entry["gram"]));m=rational_wire(entry["minimum"])
        ensure(m>0 && positive_definite(A),"form is not totally positive")
        F=form_from_trace(P,A/m);vectors=checked_vectors(P,entry["minimal_vectors"],A,m)
        V=vectors.vectors;rays=evaluation_rays(P,V)
        ensure(rank(rays)==length(P.form_basis),"form is not perfect relative to the specified space")
        index=Dict(Tuple(v)=>j for (j,v) in enumerate(V))
        push!(signed_index_maps,[index[signed_key(v)] for v in vectors.supplied])
        gens=[internal_matrix(P,matrix_wire(G)) for G in entry["stabilizer_generators_row"]]
        for g in gens;ensure(act_form(P,F,g)==F,"backend stabilizer fails exact congruence");end
        claimed=parse(BigInt,entry["stabilizer_order"])
        ensure(claimed>0,"invalid stabilizer order")
        closure=finite_matrix_closure([encode_matrix(P,g) for g in gens],P.rank*P.degree,claimed)
        ensure(haskey(closure.seen,matrix_key(-identity_matrix(QQ,P.rank*P.degree))),"stabilizer lost central -I")
        independent=finite_stabilizer(P,F)
        ensure(independent.order==claimed,"backend stabilizer is not the complete Hecke K-linear stabilizer")
        permutations=[[index[signed_key([sum(M[j,k]*v[k] for k in eachindex(v)) for j in eachindex(v)])] for v in V] for M in closure.elements]
        central=-identity_matrix(P.field,P.rank);central in gens || push!(gens,central)
        push!(actions,(;closure,permutations))
        push!(facets_all,exact_facets(P,rays))
        push!(records,PerfectFormRecord(F,A/m,QQ(1),V,rays,Any[gens...],claimed,Any[]))
    end
    for i in eachindex(records),j in 1:i-1
        ensure(isnothing(equivalence(P,records[i].form,records[j].form)),"two reported orbit representatives are K-linearly equivalent")
    end
    neighbours=NeighbourRecord[];coverage=Any[]
    for (i,entry) in enumerate(raw)
        source=records[i];vectors=source.minimal_vectors;actual_facets=facets_all[i]
        facet_lookup=Dict(Tuple(f.indices)=>j for (j,f) in enumerate(actual_facets))
        covered=Dict{Tuple,Any}()
        for (a,edge) in enumerate(entry["neighbours"])
            target=Int(edge["target"]);ensure(1<=target<=length(records),"invalid target representative")
            signed=signed_index_maps[i]
            inds=Int.(edge["facet_vector_indices_zero_based"])
            ensure(!isempty(inds) && all(0<=k<length(signed) for k in inds),"invalid signed facet incidence")
            facet=sort!(unique(signed[k+1] for k in inds));key=Tuple(facet)
            ensure(haskey(facet_lookup,key),"exported adjacency is not an OSCAR facet")
            Hraw=internal_gram(P,matrix_wire(edge["neighbour_gram"]));m=rational_wire(edge["neighbour_minimum"])
            ensure(m>0,"invalid neighbouring minimum")
            hv=checked_vectors(P,edge["neighbour_vectors"],Hraw,m);H=Hraw/m
            HF=form_from_trace(P,H)
            ensure(H!=source.trace_gram && rank(evaluation_rays(P,hv.vectors))==length(P.form_basis),"neighbour is not a distinct perfect form")
            shared=findall(v->quadratic_value(H,v)==1,vectors)
            ensure(shared==facet,"claimed facet is not the complete common minimal-vector face")
            direction=form_coordinates(P,HF-source.form);normal=actual_facets[facet_lookup[key]].normal
            pivot=findfirst(!iszero,normal);scale=direction[pivot]/normal[pivot]
            ensure(scale>0 && direction==scale*normal,"neighbour is not on the outward exact facet ray")
            g=internal_matrix(P,matrix_wire(edge["transporter_row"]))
            ensure(act_form(P,records[target].form,g)==HF,"labelled transporter has wrong direction or congruence")
            R=encode_matrix(P,g);target_index=Dict(Tuple(v)=>j for (j,v) in enumerate(records[target].minimal_vectors))
            map_vector(v)=[BigInt(ZZ(sum(R[j,k]*v[k] for k in eachindex(v)))) for j in eachindex(v)]
            ensure(Set(signed_key(map_vector(v)) for v in hv.vectors)==Set(keys(target_index)),"transporter does not map the exhaustive neighbour minimum set onto its target")
            mapping=Int[]
            for k in facet
                image=map_vector(vectors[k]);idx=target_index[signed_key(image)]
                push!(mapping,Tuple(image)==signed_key(image) ? idx : -idx)
            end
            ensure(sort!(unique(abs.(mapping))) in [f.indices for f in facets_all[target]],"transported common facet is not a target facet")
            push!(neighbours,NeighbourRecord(i,facet,target,H,g,mapping))
            for (h,permutation) in enumerate(actions[i].permutations)
                image=Tuple(sort(permutation[facet]))
                ensure(haskey(facet_lookup,image),"stabilizer does not permute the exact facets")
                if haskey(covered,image)
                    ensure(covered[image].orbit==a,"duplicate backend facet-orbit representative")
                else
                    covered[image]=(;orbit=a,stabilizer_element=h)
                end
            end
        end
        ensure(Set(keys(covered))==Set(keys(facet_lookup)),"unresolved facet orbit: enumeration is not closed")
        for (j,f) in enumerate(actual_facets)
            label=covered[Tuple(f.indices)];h=label.stabilizer_element
            M=actions[i].closure.elements[h]
            inverse_permutation=invperm(actions[i].permutations[h])
            representative_edge=entry["neighbours"][label.orbit]
            representative_facet=sort!(unique(signed_index_maps[i][k+1] for k in Int.(representative_edge["facet_vector_indices_zero_based"])))
            ensure(sort(inverse_permutation[f.indices])==representative_facet,"facet reduction permutation has the wrong direction")
            push!(source.facets,(;indices=f.indices,normal=f.normal,neighbour_orbit=label.orbit,
                to_orbit_representative=decode_matrix(P,inv(M)),
                stabilizer_word=-reverse(actions[i].closure.words[h]),
                induced_minimal_vector_permutation=inverse_permutation))
        end
        push!(coverage,Dict("representative"=>i,"facets"=>length(actual_facets),"facet_orbits"=>length(entry["neighbours"]),"unresolved_facets"=>0))
    end
    # Each mapped common facet must return to the source orbit from the target.
    for edge in neighbours
        target_facet=Tuple(sort(abs.(edge.facet_map)))
        matches=[other for other in neighbours if other.source==edge.target && any(Tuple(sort(p[other.facet]))==target_facet for p in actions[edge.target].permutations)]
        ensure(length(matches)==1 && only(matches).target==edge.source,"reverse facet identification does not close")
    end
    provenance=Dict{String,Any}("backend_commit"=>data["backend_commit"],"schema"=>SCHEMA,
        "group_constraints"=>P.group_constraints,"action_convention"=>P.action_convention,
        "minimum_and_minimal_vectors"=>"exhaustive Hecke shortest vectors; exact squared trace lengths",
        "stabilizers"=>"generated closure equals independent Hecke K-linear automorphism order; central -I retained",
        "facets"=>"exact OSCAR/polymake cone facets; every facet covered by one stabilizer orbit",
        "neighbours"=>"exact common face, outward ray, endpoint minima and full minimum-set transport; reverse incidences checked",
        "closure"=>coverage,"global_argument"=>"connectedness of the arithmetic perfect-form neighbour graph; all facet orbits exhausted under the specified K-linear GL group",
        "presentation_status"=>"not constructed; no ordered two-cell attaching loops yet")
    ArithmeticVoronoiEnumeration(P,records,neighbours,Any[],:complete,provenance)
end
