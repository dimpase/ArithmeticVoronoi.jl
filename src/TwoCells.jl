export RidgeRecord, TwoCellContext, TwoCellCertificate, trace_pairing_matrix,
       enumerate_ridges, walk_ridge, validate_two_cell, two_cell_data,
       enumerate_two_cells, write_first_two_cell, eval_matrix_word

const TWO_CELL_SCHEMA="arithmetic_voronoi_two_cells_v1"
const MatrixLetter=Tuple{Symbol,Int,Int}
qmatrix_key(A)=Tuple(fraction_string(x) for x in A)
function apply_integral(M,v)
    [BigInt(ZZ(sum(M[i,j]*v[j] for j in eachindex(v)))) for i in 1:nrows(M)]
end
inverse_word(w)=reverse([-x for x in w])
inverse_matrix_word(w)=[(k,i,-j) for (k,i,j) in reverse(w)]

function trace_pairing_matrix(P)
    basis=[form_from_trace(P,B) for B in P.form_basis]
    pairing=matrix(QQ,[P.degree==1 ? tr(A*B) : tr(tr(A*B)) for A in basis,B in basis])
    ensure(det(pairing)!=0,"degenerate trace pairing")
    pairing
end

function primitive_ray(v)
    denominator_lcm=foldl(lcm,(denominator(QQ(x)) for x in v);init=ZZ(1))
    integers=[ZZ(denominator_lcm*x) for x in v]
    divisor=foldl(gcd,integers;init=ZZ(0));ensure(divisor!=0,"zero ray")
    Tuple(BigInt(div(x,divisor)) for x in integers)
end

struct RidgeRecord
    representative::Int
    id::Int
    indices::Vector{Int}
    rays::QQMatrix
    facets::Vector{Int}
    field_rank::Int
    rank_witness::Tuple{Int,Int}
    interior_form::Any
end

"""All exact codimension-two faces, without arithmetic orbit identifications."""
function enumerate_ridges(E::ArithmeticVoronoiEnumeration)
    P=E.problem;dimension=length(P.form_basis);pairing=trace_pairing_matrix(P)
    P.rank==2 || throw(ArgumentError("presentation ridge certificates currently require rank 2; use voronoi_complex for higher-rank cells"))
    result=Vector{RidgeRecord}[]
    for (i,r) in enumerate(E.representatives)
        seen=Dict{Tuple,Vector{Int}}()
        for a in eachindex(r.facets),b in 1:a-1
            indices=intersect(r.facets[a].indices,r.facets[b].indices)
            isempty(indices) && continue
            rank(r.evaluation_rays[indices,:])==dimension-2 || continue
            seen[Tuple(indices)]=indices
        end
        records=RidgeRecord[]
        for (id,key) in enumerate(sort!(collect(keys(seen))))
            indices=seen[key]
            incident=findall(f->issubset(indices,f.indices),r.facets)
            ensure(length(incident)==2,"a ridge must have exactly two irredundant incident facets")
            field_vectors=[field_vector(P,r.minimal_vectors[j]) for j in indices]
            witness=(0,0)
            for a in eachindex(field_vectors),b in 1:a-1
                field_vectors[a][1]*field_vectors[b][2]-field_vectors[a][2]*field_vectors[b][1]!=0 && (witness=(a,b);break)
            end
            kr=witness==(0,0) ? 1 : 2
            ensure(kr==2,"ridge does not meet the positive cone: supporting vectors have K-rank one")
            interior=zero_matrix(P.field,P.rank,P.rank)
            for (v,j) in zip(field_vectors,indices)
                outer=vector_outer(P,v)
                coordinates=form_coordinates(P,outer)
                dual=[sum(pairing[a,b]*coordinates[b] for b in eachindex(coordinates)) for a in eachindex(coordinates)]
                ensure(dual==[r.evaluation_rays[j,k] for k in 1:ncols(r.evaluation_rays)],"primal/dual trace pairing mismatch")
                interior+=outer
            end
            ensure(positive_definite(trace_gram(P,interior)),"ridge positive interior witness failed")
            cone=positive_hull(r.evaluation_rays[indices,:])
            extreme=matrix(QQ,rays(cone));keys=sort!(unique(primitive_ray(v) for v in matrix_rows(extreme)))
            exact_rays=matrix(QQ,length(keys),dimension,[QQ(keys[a][b]) for a in eachindex(keys) for b in 1:dimension])
            ensure(rank(exact_rays)==dimension-2,"ridge extreme-ray dimension mismatch")
            push!(records,RidgeRecord(i,id,indices,exact_rays,incident,kr,witness,interior))
        end
        push!(result,records)
    end
    result
end

struct TwoCellContext
    enumeration::ArithmeticVoronoiEnumeration
    ridges::Vector{Vector{RidgeRecord}}
    groups::Vector{Any}
    neighbours::Vector{Vector{NeighbourRecord}}
    source_hash::String
end
function TwoCellContext(E::ArithmeticVoronoiEnumeration)
    ensure(E.completion_status==:complete && isempty(E.pending_facets),"two-cells require a complete validated enumeration")
    P=E.problem
    groups=Any[finite_matrix_closure([encode_matrix(P,g) for g in r.stabilizer],P.rank*P.degree,r.stabilizer_order) for r in E.representatives]
    data=enumeration_data(E)
    digest=bytes2hex(sha256(JSON3.write((;problem=data.problem,forms=data.forms,neighbours=data.neighbours))))
    TwoCellContext(E,enumerate_ridges(E),groups,
        [[n for n in E.neighbours if n.source==i] for i in eachindex(E.representatives)],digest)
end

function compatible_two_cell_cache(ctx::TwoCellContext,directory)
    path=joinpath(directory,"coverage.json")
    isfile(path) || return false
    ledger=JSON3.read(read(path,String),Dict{String,Any})
    get(ledger,"schema",nothing)==TWO_CELL_SCHEMA && get(ledger,"status",nothing)=="complete" &&
        get(ledger,"source_hash",nothing)==ctx.source_hash
end

function eval_matrix_word(ctx::TwoCellContext,word)
    E=ctx.enumeration;P=E.problem;value=identity_matrix(P.field,P.rank)
    for (kind,i,j) in word
        ensure(kind in (:s,:b) && j!=0,"invalid matrix-word letter")
        M=kind==:s ? E.representatives[i].stabilizer[abs(j)] : ctx.neighbours[i][abs(j)].transporter
        value*=j>0 ? M : inv(M)
    end
    value
end

function stabilizer_word(ctx,i,g)
    key=matrix_key(encode_matrix(ctx.enumeration.problem,g))
    ensure(haskey(ctx.groups[i].seen,key),"residual matrix is not in the full vertex stabilizer")
    ctx.groups[i].words[ctx.groups[i].seen[key]]
end

function facet_absolute(P,r,a,facet)
    inverse=inv(encode_matrix(P,a))
    sort!(collect(signed_key(apply_integral(inverse,r.minimal_vectors[j])) for j in r.facets[facet].indices))
end
function actual_state(ctx,i,a,incoming,outgoing)
    P=ctx.enumeration.problem;r=ctx.enumeration.representatives[i]
    (;representative=i,frame=a,gram=trace_gram(P,act_form(P,r.form,a)),incoming,outgoing,
      incoming_vectors=facet_absolute(P,r,a,incoming),outgoing_vectors=facet_absolute(P,r,a,outgoing))
end
state_key(state)=(qmatrix_key(state.gram),Tuple(state.incoming_vectors),Tuple(state.outgoing_vectors))

function step_across(ctx,state,fixed_vectors)
    E=ctx.enumeration;P=E.problem;i=state.representative;r=E.representatives[i]
    f=r.facets[state.outgoing];s=f.to_orbit_representative
    edge=ctx.neighbours[i][f.neighbour_orbit];b=edge.transporter;j=edge.target
    next_frame=b*s*state.frame
    check_group_matrix(P,encode_matrix(P,next_frame))
    H=form_from_trace(P,edge.neighbour_gram)
    ensure(act_form(P,H,s*state.frame)==act_form(P,E.representatives[j].form,next_frame),"noncommutative transporter update failed")
    target=E.representatives[j];coordinate_index=Dict(Tuple(v)=>k for (k,v) in enumerate(target.minimal_vectors))
    R=encode_matrix(P,next_frame)
    local_indices(vectors)=sort!(collect(coordinate_index[signed_key(apply_integral(R,v))] for v in vectors))
    crossed=local_indices(state.outgoing_vectors)
    incoming=findall(g->g.indices==crossed,target.facets)
    ensure(length(incoming)==1,"crossed facet is not an exact target facet")
    ridge_indices=local_indices(fixed_vectors)
    matching=findall(ridge->ridge.indices==ridge_indices,ctx.ridges[j])
    ensure(length(matching)==1,"fixed absolute ridge is not a target ridge")
    ridge=ctx.ridges[j][only(matching)]
    ensure(only(incoming) in ridge.facets,"incoming facet does not contain the fixed ridge")
    outgoing=only(setdiff(ridge.facets,incoming))
    next=actual_state(ctx,j,next_frame,only(incoming),outgoing)
    ensure(next.incoming_vectors==state.outgoing_vectors,"crossed absolute facets disagree")
    ensure(issubset(Set(signed_key(v) for v in fixed_vectors),Set(state.outgoing_vectors)),"crossing lost the fixed ridge")
    letters=vcat(MatrixLetter[(:b,i,f.neighbour_orbit)],MatrixLetter[(:s,i,k) for k in f.stabilizer_word])
    ensure(eval_matrix_word(ctx,letters)==b*s,"facet stabilizer word has the wrong direction")
    (;next,s,b,letters,source_ridge_indices=ridge_indices,facet_orbit=f.neighbour_orbit)
end

struct TwoCellCertificate
    source::Int
    ridge::Int
    initial_outgoing::Int
    steps::Vector{Any}
    residual::Any
    residual_word::Vector{Int}
    attaching_word::Vector{MatrixLetter}
    source_hash::String
end

function walk_ridge(ctx::TwoCellContext,i,ridge_id;outgoing=nothing,max_steps=1000,checkpoint=nothing)
    E=ctx.enumeration;P=E.problem;ridge=ctx.ridges[i][ridge_id];r=E.representatives[i]
    outgoing=isnothing(outgoing) ? first(ridge.facets) : outgoing
    ensure(outgoing in ridge.facets,"initial outgoing facet is not incident to the ridge")
    incoming=only(setdiff(ridge.facets,[outgoing]))
    state=actual_state(ctx,i,identity_matrix(P.field,P.rank),incoming,outgoing)
    initial=state_key(state);seen=Set{Tuple}([initial]);fixed=r.minimal_vectors[ridge.indices]
    steps=Any[];word=MatrixLetter[]
    for number in 1:max_steps
        crossing=step_across(ctx,state,fixed)
        word=[crossing.letters;word]
        ensure(eval_matrix_word(ctx,word)==crossing.next.frame,"accumulated word does not evaluate to the actual frame")
        push!(steps,(;state,s=crossing.s,b=crossing.b,facet_orbit=crossing.facet_orbit,
            letters=crossing.letters,accumulated_word=copy(word),next=crossing.next))
        state=crossing.next
        if state_key(state)==initial
            ensure(state.representative==i,"geometric closure has an inequivalent representative")
            h=state.frame;hword=stabilizer_word(ctx,i,h)
            ensure(act_form(P,r.form,h)==r.form,"residual fails exact vertex stabilization")
            correction=MatrixLetter[(:s,i,k) for k in hword]
            attaching=[word;inverse_matrix_word(correction)]
            ensure(eval_matrix_word(ctx,attaching)==identity_matrix(P.field,P.rank),"residual-corrected attaching relation is not the identity")
            certificate=TwoCellCertificate(i,ridge_id,outgoing,steps,h,hword,attaching,ctx.source_hash)
            validate_two_cell(ctx,certificate)
            return certificate
        end
        if state_key(state) in seen
            isnothing(checkpoint) || atomic_json(checkpoint,Dict("schema"=>TWO_CELL_SCHEMA,"status"=>"incomplete","reason"=>"unexpected repeated state","steps"=>number))
            throw(ArithmeticValidationError("ridge walk repeated a noninitial oriented geometric state"))
        end
        push!(seen,state_key(state))
        isnothing(checkpoint) || atomic_json(checkpoint,Dict("schema"=>TWO_CELL_SCHEMA,"status"=>"partial","source_hash"=>ctx.source_hash,"source"=>i,"ridge"=>ridge_id,"steps"=>number))
    end
    throw(ArithmeticValidationError("ridge walk reached its step limit; no completed certificate"))
end

function validate_two_cell(ctx,C::TwoCellCertificate)
    ensure(C.source_hash==ctx.source_hash,"two-cell certificate belongs to another enumeration")
    ensure(!isempty(C.steps),"empty attaching loop")
    E=ctx.enumeration;P=E.problem;r=E.representatives[C.source];ridge=ctx.ridges[C.source][C.ridge]
    fixed=r.minimal_vectors[ridge.indices]
    incoming=only(setdiff(ridge.facets,[C.initial_outgoing]))
    state=actual_state(ctx,C.source,identity_matrix(P.field,P.rank),incoming,C.initial_outgoing)
    initial=state_key(state);seen=Set{Tuple}([initial]);word=MatrixLetter[]
    for (number,step) in enumerate(C.steps)
        ensure(state_key(step.state)==state_key(state) && step.state.frame==state.frame && step.state.representative==state.representative,"recorded ridge state is inconsistent")
        ensure(step.state.incoming==state.incoming && step.state.outgoing==state.outgoing,"recorded local facet IDs are inconsistent")
        fresh=step_across(ctx,state,fixed)
        ensure(step.s==fresh.s && step.b==fresh.b && step.letters==fresh.letters,"corrupt neighbour or stabilizer label")
        ensure(step.facet_orbit==fresh.facet_orbit,"corrupt directed facet-orbit identifier")
        ensure(step.next.frame==fresh.next.frame && state_key(step.next)==state_key(fresh.next),"recorded next state is inconsistent")
        ensure(step.next.incoming==fresh.next.incoming && step.next.outgoing==fresh.next.outgoing,"recorded next local facet IDs are inconsistent")
        word=[step.letters;word]
        ensure(word==step.accumulated_word && eval_matrix_word(ctx,word)==fresh.next.frame,"corrupt accumulated frame word")
        state=fresh.next
        number<length(C.steps) && ensure(!(state_key(state) in seen),"loop contains an early geometric repetition")
        push!(seen,state_key(state))
    end
    ensure(state_key(state)==initial,"attaching loop does not close geometrically")
    ensure(state.frame==C.residual,"residual differs from the closing transporter")
    correction=MatrixLetter[(:s,C.source,k) for k in C.residual_word]
    ensure(eval_matrix_word(ctx,correction)==C.residual,"corrupt residual stabilizer correction")
    expected=[word;inverse_matrix_word(correction)]
    ensure(C.attaching_word==expected && eval_matrix_word(ctx,expected)==identity_matrix(P.field,P.rank),"invalid residual-corrected relation")
    true
end

function field_matrix_wire(P,rows)
    ensure(length(rows)==P.rank && all(length(row)==P.rank for row in rows),"wrong field-matrix shape")
    entries=Any[]
    for row in rows,coefficients in row
        ensure(length(coefficients)==P.degree,"wrong number-field coefficient count")
        push!(entries,sum(P.field(rational_wire(coefficients[k+1]))*P.generator^k for k in 0:P.degree-1))
    end
    matrix(P.field,P.rank,P.rank,entries)
end
function matrix_word_wire(raw)
    MatrixLetter[(Symbol(letter[1]),Int(letter[2]),Int(letter[3])) for letter in raw]
end
function state_wire(P,s)
    (;representative=Int(s["representative"]),frame=field_matrix_wire(P,s["frame"]),gram=matrix_wire(s["trace_gram"]),
      incoming=Int(s["incoming"]),outgoing=Int(s["outgoing"]),
      incoming_vectors=[Tuple(parse.(BigInt,v)) for v in s["incoming_vectors"]],
      outgoing_vectors=[Tuple(parse.(BigInt,v)) for v in s["outgoing_vectors"]])
end

function read_two_cell(ctx,path)
    data=JSON3.read(read(path,String),Dict{String,Any});P=ctx.enumeration.problem
    ensure(get(data,"schema",nothing)==TWO_CELL_SCHEMA && get(data,"status",nothing)=="complete","incomplete or incompatible two-cell cache")
    ensure(data["source_hash"]==ctx.source_hash,"two-cell cache has the wrong source enumeration")
    steps=Any[(;state=state_wire(P,s["state"]),s=field_matrix_wire(P,s["s"]),b=field_matrix_wire(P,s["b"]),
        facet_orbit=Int(s["facet_orbit"]),letters=matrix_word_wire(s["letters"]),
        accumulated_word=matrix_word_wire(s["accumulated_word"]),next=state_wire(P,s["next"])) for s in data["steps"]]
    C=TwoCellCertificate(Int(data["source"]),Int(data["ridge"]),Int(data["initial_outgoing"]),steps,
        field_matrix_wire(P,data["residual_matrix"]),Int.(data["residual_stabilizer_word"]),matrix_word_wire(data["attaching_word"]),data["source_hash"])
    validate_two_cell(ctx,C)
    ensure(JSON3.read(JSON3.write(two_cell_data(ctx,C)),Dict{String,Any})==data,"two-cell cache metadata or exact incidence data is inconsistent")
    C
end

function reverse_cycle_check(ctx,forward,backward)
    ensure(forward.source==backward.source && forward.ridge==backward.ridge,"reverse certificates use different based ridges")
    n=length(forward.steps);ensure(length(backward.steps)==n,"reverse walk has a different cycle length")
    corrections=Any[]
    for k in 0:n-1
        f=forward.steps[mod(-k,n)+1].state;b=backward.steps[k+1].state
        ensure(f.gram==b.gram && f.representative==b.representative,"reverse traversal has a different ordered geometric cycle")
        ensure(f.incoming_vectors==b.outgoing_vectors && f.outgoing_vectors==b.incoming_vectors,"reverse traversal does not swap incoming/outgoing absolute facets")
        correction=b.frame*inv(f.frame);word=stabilizer_word(ctx,f.representative,correction)
        ensure(eval_matrix_word(ctx,MatrixLetter[(:s,f.representative,j) for j in word])*f.frame==b.frame,"reverse frame correction failed")
        push!(corrections,(;forward_step=mod(-k,n)+1,backward_step=k+1,representative=f.representative,
            matrix=field_matrix_data(ctx.enumeration.problem,correction),stabilizer_word=word))
    end
    corrections
end

function local_ridge_orbits(ctx,i)
    P=ctx.enumeration.problem;r=ctx.enumeration.representatives[i]
    vectors=Dict(Tuple(v)=>j for (j,v) in enumerate(r.minimal_vectors))
    lookup=Dict(Tuple(ridge.indices)=>ridge.id for ridge in ctx.ridges[i])
    covered=Set{Int}();orbits=Any[]
    for ridge in ctx.ridges[i]
        ridge.id in covered && continue
        members=Dict{Int,Any}()
        for (h,M) in enumerate(ctx.groups[i].elements)
            image=Tuple(sort!(collect(vectors[signed_key(apply_integral(M,r.minimal_vectors[j]))] for j in ridge.indices)))
            ensure(haskey(lookup,image),"vertex stabilizer does not permute the exact ridges")
            id=lookup[image];push!(covered,id)
            members[id]=(;ridge=id,to_representative=field_matrix_data(P,decode_matrix(P,inv(M))),
                stabilizer_word=inverse_word(ctx.groups[i].words[h]))
        end
        push!(orbits,(;representative=ridge.id,members=[members[k] for k in sort!(collect(keys(members)))]))
    end
    ensure(length(covered)==length(ctx.ridges[i]),"local ridge reduction missed a ridge")
    orbits
end

"""Compute every local ridge in both orientations, retaining validated per-loop caches."""
function enumerate_two_cells(ctx::TwoCellContext;directory=joinpath(PACKAGE_ROOT,"results/two_cells","D$(ctx.enumeration.problem.discriminant)"),stop_after=nothing)
    mkpath(directory);cache=joinpath(directory,"loops");mkpath(cache)
    ledger=Dict{String,Any}("schema"=>TWO_CELL_SCHEMA,"status"=>"partial","source_hash"=>ctx.source_hash,
        "problem"=>problem_descriptor(ctx.enumeration.problem),"local_ridge_counts"=>length.(ctx.ridges),
        "labelled_facet_count_convention"=>"local directed facet orbits under each finite perfect-form stabilizer",
        "entries"=>[Dict("source"=>i,"ridge"=>r.id,"status"=>"pending") for (i,ridges) in enumerate(ctx.ridges) for r in ridges])
    ledger_path=joinpath(directory,"coverage.json");atomic_json(ledger_path,ledger)
    certificates=TwoCellCertificate[];cached=0;start=time_ns();cursor=0;last_checkpoint=start
    try
        for (i,ridges) in enumerate(ctx.ridges),ridge in ridges
            cursor+=1;entry=ledger["entries"][cursor];based=TwoCellCertificate[];files=String[]
            for outgoing in ridge.facets
                file=joinpath(cache,"cone$(i)_ridge$(ridge.id)_out$(outgoing).json")
                if isfile(file)
                    C=read_two_cell(ctx,file);cached+=1
                else
                    C=walk_ridge(ctx,i,ridge.id;outgoing,checkpoint=file*".partial")
                    atomic_json(file,two_cell_data(ctx,C))
                end
                push!(certificates,C);push!(based,C);push!(files,relpath(file,directory))
            end
            entry["reverse_frame_corrections"]=reverse_cycle_check(ctx,based...)
            entry["status"]="complete";entry["certificates"]=files
            entry["length"]=length(first(based).steps)
            entry["identification_witness"]="identity: this local ridge is processed directly in its stored frame"
            ledger["cache_hits"]=cached;ledger["completed_oriented_loops"]=length(certificates)
            # Each loop already has an atomic certificate file, from which a
            # resumed run reconstructs coverage. Rewriting the growing ledger
            # after every ridge would make checkpoint work quadratic in the
            # number of ridges. Refresh live coverage about once per second;
            # completion, requested stops and errors always save it below.
            if cursor==1 || time_ns()-last_checkpoint>=1_000_000_000
                atomic_json(ledger_path,ledger);last_checkpoint=time_ns()
            end
            if cursor==1 && i==1
                write_first_two_cell(ctx,first(based);directory)
            end
            println("D",ctx.enumeration.problem.discriminant," ridge ",cursor,"/",length(ledger["entries"]),": cycle length ",entry["length"]);flush(stdout)
            if !isnothing(stop_after) && cursor>=stop_after
                ledger["stop_reason"]="requested useful-batch limit; remaining entries stay pending"
                atomic_json(ledger_path,ledger)
                return (;certificates,ledger)
            end
        end
        ledger["local_ridge_orbits"]=[local_ridge_orbits(ctx,i) for i in eachindex(ctx.ridges)]
        ledger["locally_reduced_ridge_counts"]=length.(ledger["local_ridge_orbits"])
        ledger["completed_oriented_loops"]=length(certificates)
        ledger["loop_lengths"]=[length(C.steps) for C in certificates]
        ledger["global_two_cell_orbits"]=nothing
        ledger["global_reduction"]="not performed; all local ridges and both orientations are retained as redundant attaching loops"
        ledger["cache_hits"]=cached;ledger["elapsed_seconds"]=(time_ns()-start)/1e9
        ledger["status"]="complete";atomic_json(ledger_path,ledger)
    catch err
        ledger["status"]="incomplete";ledger["error"]=sprint(showerror,err);atomic_json(ledger_path,ledger)
        rethrow()
    end
    (;certificates,ledger)
end

matrix_word_data(word)=[[string(k),i,j] for (k,i,j) in word]
function state_data(P,s)
    (;representative=s.representative,frame=field_matrix_data(P,s.frame),trace_gram=matrix_strings(s.gram),
      incoming=s.incoming,outgoing=s.outgoing,incoming_vectors=[[string(x) for x in v] for v in s.incoming_vectors],
      outgoing_vectors=[[string(x) for x in v] for v in s.outgoing_vectors])
end
function ridge_data(ctx,ridge)
    P=ctx.enumeration.problem;r=ctx.enumeration.representatives[ridge.representative]
    (;representative=ridge.representative,id=ridge.id,dimension=length(P.form_basis)-2,
      minimal_vector_indices=ridge.indices,extreme_rays=matrix_strings(ridge.rays),incident_facets=ridge.facets,
      supporting_vectors=[[string(x) for x in v] for v in r.minimal_vectors[ridge.indices]],
      K_rank=ridge.field_rank,rank_witness=collect(ridge.rank_witness),
      positive_interior_form=field_matrix_data(P,ridge.interior_form))
end
function two_cell_data(ctx,C)
    P=ctx.enumeration.problem
    (;schema=TWO_CELL_SCHEMA,status="complete",source_hash=C.source_hash,problem=problem_descriptor(P),
      group=P.group_constraints,action=P.action_convention,
      word_convention="written matrix-product order; steps prepend b*s; s indexes vertex stabilizers, b indexes local directed facet orbits; negative index means inverse",
      pairing_matrix=matrix_strings(trace_pairing_matrix(P)),source=C.source,ridge=C.ridge,
      ridge_data=ridge_data(ctx,ctx.ridges[C.source][C.ridge]),initial_outgoing=C.initial_outgoing,
      starting_form=field_matrix_data(P,ctx.enumeration.representatives[C.source].form),
      steps=[(;state=state_data(P,s.state),s=field_matrix_data(P,s.s),b=field_matrix_data(P,s.b),
        facet_orbit=s.facet_orbit,letters=matrix_word_data(s.letters),accumulated_word=matrix_word_data(s.accumulated_word),
        next=state_data(P,s.next)) for s in C.steps],
      residual_matrix=field_matrix_data(P,C.residual),residual_stabilizer_word=C.residual_word,
      attaching_word=matrix_word_data(C.attaching_word),closed_geometrically=true,matrix_relation_verified=true)
end

function write_first_two_cell(ctx,C;directory=joinpath(PACKAGE_ROOT,"results/two_cells","D$(ctx.enumeration.problem.discriminant)"))
    mkpath(directory);atomic_json(joinpath(directory,"first_two_cell.json"),two_cell_data(ctx,C))
    P=ctx.enumeration.problem;r=ctx.ridges[C.source][C.ridge]
    open(joinpath(directory,"FIRST_TWO_CELL.txt"),"w") do io
        println(io,"First exact two-cell, field discriminant ",P.discriminant)
        println(io,"Group: ",P.group_constraints,". Scalar units and -I retained.")
        println(io,"Source form ",C.source,"; ridge ",C.ridge,"; incident facets ",r.facets,"; K-rank ",r.field_rank)
        println(io,"Dual/primal trace pairing matrix: ",matrix_strings(trace_pairing_matrix(P)))
        println(io,"Starting form coefficients: ",field_matrix_data(P,ctx.enumeration.representatives[C.source].form))
        println(io,"Ridge extreme rays: ",matrix_strings(r.rays))
        println(io,"Ridge supporting vectors: ",ctx.enumeration.representatives[C.source].minimal_vectors[r.indices])
        println(io,"Rank witness positions: ",r.rank_witness,"; positive interior form: ",field_matrix_data(P,r.interior_form))
        for (k,step) in enumerate(C.steps)
            println(io,"Step ",k,": representative ",step.state.representative,"; incoming/outgoing facets ",step.state.incoming,"/",step.state.outgoing)
            println(io,"  a = ",field_matrix_data(P,step.state.frame))
            println(io,"  s = ",field_matrix_data(P,step.s),"; b = ",field_matrix_data(P,step.b))
            println(io,"  next a = b*s*a = ",field_matrix_data(P,step.next.frame))
            println(io,"  accumulated word = ",matrix_word_data(step.accumulated_word))
        end
        println(io,"Residual h = ",field_matrix_data(P,C.residual),"; finite stabilizer word ",C.residual_word)
        println(io,"Corrected attaching word: ",matrix_word_data(C.attaching_word))
        println(io,"Verified: actual oriented state closes; the word evaluates to h; correction evaluates to h^-1; corrected relation is I.")
    end
end
