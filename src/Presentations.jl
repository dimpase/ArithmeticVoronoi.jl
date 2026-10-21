export PresentationCertificate, build_gl_presentation, presentation_data,
       verify_presentation, matrix_to_word

const PRESENTATION_SCHEMA="arithmetic_voronoi_presentation_v1"
function freely_reduce(word)
    result=Int[]
    for letter in word
        if !isempty(result) && last(result)==-letter;pop!(result)
        else;push!(result,letter)
        end
    end
    result
end
function evaluate_word(matrices,word;identity=nothing)
    value=isnothing(identity) ? identity_matrix(base_ring(first(matrices)),nrows(first(matrices))) : identity
    for letter in word;value*=letter>0 ? matrices[letter] : inv(matrices[-letter]);end
    value
end
vertex_elements(ctx,i)=[decode_matrix(ctx.enumeration.problem,M) for M in ctx.groups[i].elements]

"""Exact global unoriented edges and every local directed incidence orbit."""
function graph_edges(ctx)
    E=ctx.enumeration;P=E.problem;I=identity_matrix(P.field,P.rank)
    vertices=[vertex_elements(ctx,i) for i in eachindex(E.representatives)]
    edges=Any[];darts=Any[]
    for (i,list) in enumerate(ctx.neighbours),(k,n) in enumerate(list)
        F=E.representatives[i].form;H=form_from_trace(P,n.neighbour_gram)
        found=nothing
        for (q,edge) in enumerate(edges)
            A=E.representatives[edge.source].form;B=edge.neighbour
            for endpoint in 1:2
                (endpoint==1 ? edge.source : edge.target)==i || continue
                for u in vertices[i]
                    p=endpoint==1 ? u : inv(edge.transporter)*u
                    left,right=endpoint==1 ? (A,B) : (B,A)
                    if act_form(P,left,p)==F && act_form(P,right,p)==H
                        found=(q,p,endpoint);break
                    end
                end
                isnothing(found) || break
            end
            isnothing(found) || break
        end
        if isnothing(found)
            q=length(edges)+1
            push!(edges,(;id=q,source=i,target=n.target,transporter=n.transporter,neighbour=H,
                midpoint=(F+H)/2,darts=Int[]))
            found=(q,I,1)
        end
        q,p,endpoint=found;id=length(darts)+1
        push!(edges[q].darts,id)
        push!(darts,(;id,source=i,local_orbit=k,target=n.target,transporter=n.transporter,
            neighbour=H,edge=q,edge_to_local=p,source_endpoint=endpoint))
    end
    for (q,edge) in enumerate(edges)
        F=E.representatives[edge.source].form;H=edge.neighbour
        oriented=[u for u in vertices[edge.source] if act_form(P,H,u)==H]
        reversals=edge.source==edge.target ? [u*edge.transporter for u in vertices[edge.source] if act_form(P,H,u*edge.transporter)==F] : Any[]
        full=[oriented;reversals]
        ensure(length(full)==length(oriented)*(isempty(reversals) ? 1 : 2),"setwise/oriented edge stabilizer index is not one or two")
        ensure(length(edge.darts)==(isempty(reversals) ? 2 : 1),"local directed incidence orbits do not match global edge/inversion classification")
        keys=Set(matrix_key(encode_matrix(P,g)) for g in oriented)
        for g in full
            check_group_matrix(P,encode_matrix(P,g))
            ensure(act_form(P,edge.midpoint,g)==edge.midpoint,"edge stabilizer does not fix its midpoint")
        end
        for g in reversals
            ensure(act_form(P,F,g)==H && act_form(P,H,g)==F,"endpoint reversal failed")
            ensure(matrix_key(encode_matrix(P,g*g)) in keys,"endpoint reversal square is outside the oriented edge stabilizer")
        end
        edges[q]=merge(edge,(;oriented_stabilizer=oriented,setwise_stabilizer=full,reversals))
    end
    for d in darts
        E0=[u for u in vertices[d.source] if act_form(P,d.neighbour,u)==d.neighbour]
        canonical=edges[d.edge]
        Hkeys=Set(matrix_key(encode_matrix(P,g)) for g in canonical.setwise_stabilizer)
        for u in E0
            ensure(matrix_key(encode_matrix(P,d.edge_to_local*u*inv(d.edge_to_local))) in Hkeys,"half-edge inclusion into midpoint stabilizer failed")
            v=d.transporter*u*inv(d.transporter)
            ensure(act_form(P,E.representatives[d.target].form,v)==E.representatives[d.target].form,"original oriented-edge inclusion failed")
        end
        darts[d.id]=merge(d,(;oriented_stabilizer=E0))
    end
    (;edges,darts)
end

function coherent_lifts(ctx,graph)
    P=ctx.enumeration.problem;n=length(ctx.enumeration.representatives);total=n+length(graph.edges)
    lifts=Any[nothing for _ in 1:total];lifts[1]=identity_matrix(P.field,P.rank)
    tree=Set{Int}();changed=true
    while changed
        changed=false
        for d in graph.darts
            source=d.source;target=n+d.edge;p=d.edge_to_local
            if !isnothing(lifts[source]) && isnothing(lifts[target])
                lifts[target]=p*lifts[source];push!(tree,d.id);changed=true
            elseif isnothing(lifts[source]) && !isnothing(lifts[target])
                lifts[source]=inv(p)*lifts[target];push!(tree,d.id);changed=true
            end
        end
    end
    ensure(all(!isnothing,lifts) && length(tree)==total-1,"quotient incidence graph is disconnected or spanning-tree lifts are incomplete")
    labels=Any[inv(lifts[n+d.edge])*d.edge_to_local*lifts[d.source] for d in graph.darts]
    for id in tree;ensure(labels[id]==identity_matrix(P.field,P.rank),"tree label is nontrivial after choosing coherent lifts");end
    (;lifts,tree,labels)
end

function finite_group_record(P,all_elements;preferred=nothing)
    I=identity_matrix(P.field,P.rank);identity_key=matrix_key(encode_matrix(P,I))
    elements=Any[I];seen=Set([identity_key])
    for g in all_elements
        key=matrix_key(encode_matrix(P,g))
        if !(key in seen);push!(elements,g);push!(seen,key);end
    end
    keys=Dict(matrix_key(encode_matrix(P,g))=>i for (i,g) in enumerate(elements))
    function closure(gens)
        values=Any[I];words=[Int[]];known=Dict(identity_key=>1);cursor=1
        while cursor<=length(values)
            for (k,g) in enumerate(gens)
                value=values[cursor]*g;key=matrix_key(encode_matrix(P,value))
                ensure(haskey(keys,key),"finite stabilizer multiplication left the claimed complete group")
                if !haskey(known,key)
                    push!(values,value);push!(words,[words[cursor];k]);known[key]=length(values)
                end
            end
            cursor+=1
        end
        (;values,words,known)
    end
    generators=isnothing(preferred) ? Any[-I] : Any[preferred...]
    -I in generators || push!(generators,-I)
    filter!(g->g!=I,generators);unique!(generators)
    enumeration=closure(generators)
    for g in elements
        key=matrix_key(encode_matrix(P,g))
        if !haskey(enumeration.known,key)
            push!(generators,g);enumeration=closure(generators)
        end
    end
    ensure(length(enumeration.values)==length(elements),"finite group generators are incomplete")
    words=[enumeration.words[enumeration.known[matrix_key(encode_matrix(P,g))]] for g in elements]
    table=[keys[matrix_key(encode_matrix(P,g*h))] for g in elements,h in elements]
    gen_indices=[keys[matrix_key(encode_matrix(P,g))] for g in generators]
    (;elements,generators,words,keys,multiplication_table=table,generator_indices=gen_indices)
end

struct PresentationCertificate
    context::TwoCellContext
    graph::Any
    coherent::Any
    vertex_groups::Vector{Any}
    generator_matrices::Vector{Any}
    generator_provenance::Vector{Any}
    relators::Vector{Vector{Int}}
    relator_provenance::Vector{Any}
    vertex_generator_ids::Vector{Vector{Int}}
    half_edge_words::Vector{Vector{Int}}
    transporter_words::Vector{Vector{Int}}
    loop_count::Int
end

function finite_element_word(C::PresentationCertificate,vertex,g)
    P=C.context.enumeration.problem;group=C.vertex_groups[vertex]
    key=matrix_key(encode_matrix(P,g));ensure(haskey(group.keys,key),"matrix is outside the claimed finite vertex group")
    local_word=group.words[group.keys[key]]
    [j>0 ? C.vertex_generator_ids[vertex][j] : -C.vertex_generator_ids[vertex][-j] for j in local_word]
end

function build_gl_presentation(ctx::TwoCellContext,cells;coverage)
    ensure(coverage["status"]=="complete" && all(e["status"]=="complete" for e in coverage["entries"]),"cannot certify a presentation with incomplete ridge coverage")
    ensure(coverage["source_hash"]==ctx.source_hash,"ridge coverage belongs to a different source enumeration")
    expected=Set((i,r.id) for (i,list) in enumerate(ctx.ridges) for r in list)
    ensure(Set((C.source,C.ridge) for C in cells)==expected,"attaching loops do not cover every local ridge")
    for C in cells;validate_two_cell(ctx,C);end
    E=ctx.enumeration;P=E.problem;n=length(E.representatives)
    graph=graph_edges(ctx);coherent=coherent_lifts(ctx,graph)
    groups=Any[]
    for (i,r) in enumerate(E.representatives)
        a=coherent.lifts[i]
        push!(groups,finite_group_record(P,[inv(a)*u*a for u in vertex_elements(ctx,i)];preferred=[inv(a)*u*a for u in r.stabilizer]))
    end
    for edge in graph.edges
        a=coherent.lifts[n+edge.id]
        push!(groups,finite_group_record(P,[inv(a)*u*a for u in edge.setwise_stabilizer]))
    end
    matrices=Any[];origins=Any[];vertex_ids=Vector{Int}[]
    for (vertex,group) in enumerate(groups)
        ids=Int[]
        for (j,g) in enumerate(group.generators)
            push!(matrices,g);push!(ids,length(matrices));push!(origins,(;kind="finite_vertex_stabilizer",vertex,local_generator=j))
        end
        push!(vertex_ids,ids)
    end
    half_words=Vector{Int}[]
    for d in graph.darts
        if d.id in coherent.tree;push!(half_words,Int[])
        else
            push!(matrices,coherent.labels[d.id]);push!(half_words,[length(matrices)])
            push!(origins,(;kind="non_tree_half_edge",dart=d.id))
        end
    end
    relators=Vector{Int}[];provenance=Any[]
    C=PresentationCertificate(ctx,graph,coherent,groups,matrices,origins,relators,provenance,
        vertex_ids,half_words,Vector{Int}[],length(cells))
    function add_relation(word,origin)
        reduced=freely_reduce(word)
        ensure(evaluate_word(matrices,reduced)==identity_matrix(P.field,P.rank),"proposed presentation relator fails exact matrix evaluation")
        isempty(reduced) && return
        push!(relators,reduced);push!(provenance,origin)
    end
    for (vertex,group) in enumerate(groups),(x,element) in enumerate(group.elements),(j,g) in enumerate(group.generators)
        y=group.multiplication_table[x,group.generator_indices[j]]
        w=[vertex_ids[vertex][k] for k in group.words[x]]
        v=[vertex_ids[vertex][k] for k in group.words[y]]
        add_relation([w;vertex_ids[vertex][j];inverse_word(v)],(;kind="finite_group_cayley_relation",vertex,element=x,generator=j,product=y))
    end
    for d in graph.darts,u in d.oriented_stabilizer
        a=coherent.lifts[d.source];middle=n+d.edge;z=coherent.lifts[middle]
        source_element=inv(a)*u*a
        target_element=inv(z)*d.edge_to_local*u*inv(d.edge_to_local)*z
        w=finite_element_word(C,d.source,source_element);v=finite_element_word(C,middle,target_element)
        t=half_words[d.id]
        add_relation([t;w;inverse_word(t);inverse_word(v)],(;kind="half_edge_stabilizer_inclusion",dart=d.id,
            source_element=field_matrix_data(P,source_element),target_element=field_matrix_data(P,target_element)))
    end
    for d in graph.darts
        j=d.target;b=d.transporter;p=d.edge_to_local
        target_indices=sort!(abs.(ctx.neighbours[d.source][d.local_orbit].facet_map))
        target_facet=only(findall(f->f.indices==target_indices,E.representatives[j].facets))
        f=E.representatives[j].facets[target_facet];s=f.to_orbit_representative
        reverse_id=only(x.id for x in graph.darts if x.source==j && x.local_orbit==f.neighbour_orbit)
        reverse=graph.darts[reverse_id]
        ensure(reverse.edge==d.edge,"reversed edge has a different global orbit")
        z=reverse.edge_to_local*s*b*inv(p)
        middle=n+d.edge;lmid=coherent.lifts[middle];lj=coherent.lifts[j]
        zword=finite_element_word(C,middle,inv(lmid)*z*lmid)
        sword=finite_element_word(C,j,inv(lj)*s*lj)
        word=freely_reduce([inverse_word(sword);inverse_word(half_words[reverse_id]);zword;half_words[d.id]])
        expected_matrix=inv(lj)*b*coherent.lifts[d.source]
        ensure(evaluate_word(matrices,word)==expected_matrix,"rewritten original transporter has the wrong coherent lift or multiplication order")
        push!(C.transporter_words,word)
    end
    lookup=Dict((d.source,d.local_orbit)=>d.id for d in graph.darts)
    for (id,cell) in enumerate(cells)
        word=Int[]
        for step in cell.steps
            i=step.state.representative;l=coherent.lifts[i]
            s=finite_element_word(C,i,inv(l)*step.s*l)
            d=lookup[(i,step.facet_orbit)]
            word=freely_reduce([C.transporter_words[d];s;word])
        end
        i=cell.source;l=coherent.lifts[i]
        residual=finite_element_word(C,i,inv(l)*cell.residual*l)
        add_relation([word;inverse_word(residual)],(;kind="ordered_two_cell",cell=id,source=i,ridge=cell.ridge,
            initial_outgoing=cell.initial_outgoing,length=length(cell.steps),residual=field_matrix_data(P,cell.residual)))
    end
    verify_presentation(C)
    C
end

function verify_presentation(C::PresentationCertificate)
    P=C.context.enumeration.problem;I=identity_matrix(P.field,P.rank)
    for (i,g) in enumerate(C.generator_matrices)
        check_group_matrix(P,encode_matrix(P,g));origin=C.generator_provenance[i]
        expected=origin.kind=="finite_vertex_stabilizer" ? C.vertex_groups[origin.vertex].generators[origin.local_generator] : C.coherent.labels[origin.dart]
        ensure(g==expected,"generator matrix differs from its finite-group or coherent-edge provenance")
    end
    for group in C.vertex_groups
        ensure(-I in group.generators,"finite vertex group lost the central -I generator")
    end
    for word in C.relators;ensure(evaluate_word(C.generator_matrices,word)==I,"presentation matrix relator failed");end
    ensure(C.loop_count>=sum(length, C.context.ridges),"insufficient attaching-loop coverage")
    true
end

function presentation_data(C::PresentationCertificate)
    ctx=C.context;P=ctx.enumeration.problem;n=length(ctx.enumeration.representatives)
    (;schema=PRESENTATION_SCHEMA,status="certified_GL",group=P.group_constraints,source_hash=ctx.source_hash,
      action=P.action_convention,method="Brown presentation on the edge-subdivided well-rounded complex; complete setwise edge stabilizers at new vertices",
      generators=[(;id=i,matrix=field_matrix_data(P,g),provenance=C.generator_provenance[i]) for (i,g) in enumerate(C.generator_matrices)],
      relators=[(;word=w,provenance=C.relator_provenance[i]) for (i,w) in enumerate(C.relators)],
      vertex_groups=[(;vertex=i,order=length(g.elements),elements=field_matrix_data.(Ref(P),g.elements),
          generator_ids=C.vertex_generator_ids[i],element_words=g.words,multiplication_table=[collect(row) for row in eachrow(g.multiplication_table)]) for (i,g) in enumerate(C.vertex_groups)],
      global_original_edges=[(;id=e.id,source=e.source,target=e.target,transporter=field_matrix_data(P,e.transporter),
          midpoint=field_matrix_data(P,e.midpoint),local_directed_orbits=e.darts,oriented_stabilizer=field_matrix_data.(Ref(P),e.oriented_stabilizer),
          setwise_stabilizer=field_matrix_data.(Ref(P),e.setwise_stabilizer),endpoint_reversals=field_matrix_data.(Ref(P),e.reversals),
          reversal_squares=field_matrix_data.(Ref(P),[g*g for g in e.reversals])) for e in C.graph.edges],
      half_edges=[(;id=d.id,source=d.source,target=n+d.edge,local_orbit=d.local_orbit,global_original_edge=d.edge,
          edge_to_local=field_matrix_data(P,d.edge_to_local),source_endpoint=d.source_endpoint,
          stabilizer=field_matrix_data.(Ref(P),d.oriented_stabilizer),coherent_transporter=field_matrix_data(P,C.coherent.labels[d.id]),
          word=C.half_edge_words[d.id],original_transporter_word=C.transporter_words[d.id]) for d in C.graph.darts],
      coherent_vertex_lifts=field_matrix_data.(Ref(P),C.coherent.lifts),tree_half_edges=sort!(collect(C.coherent.tree)),
      attaching_loop_count=C.loop_count,all_relators_verified=true,simplification="raw presentation retained; no generator elimination yet")
end

"""Constructive word by exact Voronoi-cone descent toward an interior target.

The fixed positive target is g^-1*B*g^-t, where B is a positive sum of every
root-cone evaluation generator. Crossing a separating facet strictly decreases
Tr(F_current*target). Values have a fixed denominator, so descent terminates.
"""
function matrix_to_word(C::PresentationCertificate,g;max_cones=10000,timeout=120.0,checkpoint=nothing)
    ctx=C.context;E=ctx.enumeration;P=E.problem;check_group_matrix(P,encode_matrix(P,g))
    ensure(C.coherent.lifts[1]==identity_matrix(P.field,P.rank),"word routine requires the root lift to be identity")
    base=zero_matrix(P.field,P.rank,P.rank)
    for z in E.representatives[1].minimal_vectors
        v=field_vector(P,z);base+=vector_outer(P,v)
    end
    fixed=inv(g)*base*form_adjoint(P,inv(g))
    ensure(positive_definite(trace_gram(P,fixed)),"constructive reduction target is not positive definite")
    targetkey=qmatrix_key(trace_gram(P,act_form(P,E.representatives[1].form,g)))
    pairing=trace_pairing_matrix(P);i=1;a=identity_matrix(P.field,P.rank);word=Int[];start=time_ns();steps=Any[]
    objective(F)=sum(form_coordinates(P,F)[k]*pairing[k,l]*form_coordinates(P,fixed)[l] for k in 1:nrows(pairing),l in 1:ncols(pairing))
    lookup=Dict((d.source,d.local_orbit)=>d.id for d in C.graph.darts)
    denominator_forms=foldl(lcm,(denominator(QQ(coeff(F.form[x,y],k))) for F in E.representatives for x in 1:P.rank for y in 1:P.rank for k in 0:P.degree-1);init=ZZ(1))
    denominator_target=foldl(lcm,(denominator(QQ(coeff(fixed[x,y],k))) for x in 1:P.rank for y in 1:P.rank for k in 0:P.degree-1);init=ZZ(1))
    objective_denominator=denominator_forms*denominator_target
    function incomplete(reason)
        isnothing(checkpoint) || atomic_json(checkpoint,(;schema="av_word_path_v1",status="incomplete",reason,
            source_hash=ctx.source_hash,input_matrix=field_matrix_data(P,g),steps,word,
            current_representative=i,current_frame=field_matrix_data(P,a),objective_denominator=string(objective_denominator)))
        throw(ArithmeticValidationError(reason))
    end
    for number in 0:max_cones-1
        local_target=a*fixed*form_adjoint(P,a);coordinates=form_coordinates(P,local_target)
        dual=[sum(pairing[k,l]*coordinates[l] for l in eachindex(coordinates)) for k in eachindex(coordinates)]
        violations=[sum(f.normal[k]*dual[k] for k in eachindex(dual)) for f in E.representatives[i].facets]
        if all(x>=0 for x in violations)
            ensure(qmatrix_key(trace_gram(P,act_form(P,E.representatives[i].form,a)))==targetkey,"interior cone reduction ended at the wrong actual perfect form")
            ensure(i==1,"target perfect form has the wrong arithmetic orbit")
            h=a*inv(g);hword=finite_element_word(C,1,h)
            result=freely_reduce([inverse_word(hword);word])
            ensure(evaluate_word(C.generator_matrices,result)==g,"constructive word does not evaluate to its input matrix")
            return (;word=result,method="exact separating-facet cone descent",visited=number+1,discovered=number+1,
                steps,residual_stabilizer_word=hword,objective_denominator=string(objective_denominator),seconds=(time_ns()-start)/1e9)
        end
        (time_ns()-start)/1e9<timeout || incomplete("constructive reduction reached its explicit time limit; word remains unresolved")
        facet=findfirst(x->x<0,violations);f=E.representatives[i].facets[facet]
        d=lookup[(i,f.neighbour_orbit)];edge=ctx.neighbours[i][f.neighbour_orbit];s=f.to_orbit_representative
        next=edge.transporter*s*a;j=edge.target
        before=objective(act_form(P,E.representatives[i].form,a));after=objective(act_form(P,E.representatives[j].form,next))
        ensure(after<before,"separating-facet crossing did not strictly decrease the exact trace objective")
        ensure(denominator(before*objective_denominator)==1 && denominator(after*objective_denominator)==1,"trace objectives left their fixed denominator lattice")
        l=C.coherent.lifts[i];sword=finite_element_word(C,i,inv(l)*s*l)
        word=freely_reduce([C.transporter_words[d];sword;word])
        ensure(evaluate_word(C.generator_matrices,word)==inv(C.coherent.lifts[j])*next,"constructive descent frame word failed")
        push!(steps,(;source=i,target=j,facet,before=fraction_string(before),after=fraction_string(after),
            source_frame=field_matrix_data(P,a),source_correction=field_matrix_data(P,s),
            neighbour_transporter=field_matrix_data(P,edge.transporter),actual_frame=field_matrix_data(P,next),word=copy(word)))
        i=j;a=next
    end
    incomplete("constructive cone descent reached its explicit step limit; word remains unresolved")
end
