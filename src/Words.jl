export WordReductionError, verify_word_certificate

struct WordReductionError <: Exception
    kind::Symbol
    message::String
    checkpoint::Union{Nothing,String}
end
Base.showerror(io::IO,e::WordReductionError)=print(io,e.kind,": ",e.message)

function word_input(R,g,target)
    P=R.enumeration.problem
    target in (:GL,:SL,:PSL) || throw(WordReductionError(:invalid_input,"target must be GL, SL or PSL",nothing))
    R isa GLGroupResult && target!=:GL && throw(WordReductionError(:invalid_input,"this result contains only a GL presentation",nothing))
    try
        value=size(g)==(P.rank*P.degree,P.rank*P.degree) ? decode_matrix(P,g) : matrix(P.field,g)
        check_group_matrix(P,encode_matrix(P,value))
        target==:GL || det(value)==1 || throw(GroupConstraintError("SL/PSL inputs must be determinant-one lifts"))
        value
    catch err
        throw(WordReductionError(:invalid_input,sprint(showerror,err),nothing))
    end
end

function word_alphabet(R,target,simplified)
    simplified && !haskey(R.simplified,target) && throw(WordReductionError(:unavailable_alphabet,"simplification for the requested target is incomplete",nothing))
    raw=target==:GL ? R.gl : (target==:SL ? R.sl : R.psl)
    object=simplified ? R.simplified[target] : raw
    payload=(;source_hash=R.context.source_hash,target=string(target),simplified,
        matrices=field_matrix_data.(Ref(R.enumeration.problem),object.generator_matrices),relators=object.relators)
    (;matrices=object.generator_matrices,fingerprint=bytes2hex(sha256(JSON3.write(payload))))
end

function translated_word(R,raw_word,target,simplified)
    word=raw_word
    if target!=:GL
        rewrite=schreier_rewrite(R.sl,raw_word)
        ensure(rewrite.coset==0,"determinant-one word has a nontrivial determinant quotient coset")
        word=rewrite.word
    end
    simplified && (word=substitute_word(word,R.simplified[target].forward_generator_words))
    freely_reduce(word)
end

function reduction_target(C,g)
    E=C.context.enumeration;P=E.problem;B=zero_matrix(P.field,P.rank,P.rank)
    for z in E.representatives[1].minimal_vectors
        v=field_vector(P,z);B+=vector_outer(P,v)
    end
    inv(g)*B*form_adjoint(P,inv(g))
end

"""Replay a supplied path without running a search or trusting its matrix word."""
function replay_word_path(C,g,data)
    E=C.context.enumeration;P=E.problem;ctx=C.context;pairing=trace_pairing_matrix(P)
    Q=reduction_target(C,g);i=1;a=identity_matrix(P.field,P.rank);word=Int[]
    lookup=Dict((d.source,d.local_orbit)=>d.id for d in C.graph.darts)
    objective(F)=sum(form_coordinates(P,F)[k]*pairing[k,l]*form_coordinates(P,Q)[l] for k in 1:nrows(pairing),l in 1:ncols(pairing))
    nf=foldl(lcm,(matrix_coefficient_denominator(P,r.form) for r in E.representatives);init=ZZ(1))
    nq=matrix_coefficient_denominator(P,Q)
    ensure(parse(BigInt,data["objective_denominator"])==nf*nq,"corrupt objective denominator bound")
    for step in data["steps"]
        ensure(step["source"]==i && field_matrix_wire(P,step["source_frame"])==a,"word path has a corrupt source frame")
        facet=Int(step["facet"]);f=E.representatives[i].facets[facet]
        localQ=a*Q*form_adjoint(P,a);coordinates=form_coordinates(P,localQ)
        dual=[sum(pairing[k,l]*coordinates[l] for l in eachindex(coordinates)) for k in eachindex(coordinates)]
        ensure(sum(f.normal[k]*dual[k] for k in eachindex(dual))<0,"word path did not cross a separating facet")
        edge=ctx.neighbours[i][f.neighbour_orbit];s=f.to_orbit_representative;j=edge.target
        ensure(field_matrix_wire(P,step["source_correction"])==s && field_matrix_wire(P,step["neighbour_transporter"])==edge.transporter,"corrupt word crossing labels")
        next=edge.transporter*s*a
        ensure(step["target"]==j && field_matrix_wire(P,step["actual_frame"])==next,"word path has a corrupt next frame")
        before=objective(act_form(P,E.representatives[i].form,a));after=objective(act_form(P,E.representatives[j].form,next))
        ensure(before==rational_wire(step["before"]) && after==rational_wire(step["after"]) && after<before,"corrupt or nondecreasing trace objectives")
        ensure(denominator(before*nf*nq)==1 && denominator(after*nf*nq)==1,"word path left its fixed rational denominator lattice")
        l=C.coherent.lifts[i];sword=finite_element_word(C,i,inv(l)*s*l)
        word=freely_reduce([C.transporter_words[lookup[(i,f.neighbour_orbit)]];sword;word])
        ensure(Int.(step["word"])==word,"corrupt word path accumulated word")
        ensure(evaluate_word(C.generator_matrices,word)==inv(C.coherent.lifts[j])*next,"word path matrix frame mismatch")
        i=j;a=next
    end
    ensure(i==1 && act_form(P,E.representatives[1].form,a)==act_form(P,E.representatives[1].form,g),"word path ends at the wrong target cone")
    h=a*inv(g);hword=Int.(data["residual_stabilizer_word"])
    ensure(all(abs(k) in C.vertex_generator_ids[1] for k in hword),"word path correction is not in the root finite-stabilizer alphabet")
    ensure(haskey(C.vertex_groups[1].keys,matrix_key(encode_matrix(P,h))),"word path residual is outside the root stabilizer")
    ensure(evaluate_word(C.generator_matrices,hword)==h,"word path residual stabilizer correction is corrupt")
    result=freely_reduce([inverse_word(hword);word])
    ensure(result==Int.(data["word"]) && evaluate_word(C.generator_matrices,result)==g,"word path final matrix word is corrupt")
    true
end

function verify_word_certificate(R,g,data)
    target=Symbol(data["target"]);simplified=Bool(data["simplified"]);P=R.enumeration.problem
    g=word_input(R,g,target);alphabet=word_alphabet(R,target,simplified)
    ensure(data["schema"] in ("av_constructive_word_v1","av_constructive_word_v2") && data["status"]=="complete","incomplete or incompatible word certificate")
    ensure(data["source_hash"]==R.context.source_hash && data["alphabet_fingerprint"]==alphabet.fingerprint,"word certificate has a different source or alphabet")
    ensure(field_matrix_wire(P,data["input_matrix"])==g,"word certificate belongs to another input matrix")
    replay_word_path(R.gl,g,data["raw_GL_path"])
    word=translated_word(R,Int.(data["raw_GL_path"]["word"]),target,simplified)
    ensure(word==Int.(data["word"]),"word certificate alphabet translation is corrupt")
    I=identity_matrix(P.field,P.rank);value=evaluate_word(alphabet.matrices,word)
    scalar=target==:PSL ? projective_scalar(value,g) : (value==g ? one(P.field) : nothing)
    ensure(!isnothing(scalar),"word certificate evaluates to the wrong matrix or scalar centre element")
    if data["schema"]=="av_constructive_word_v2"
        ensure(data["SL_lift_scalar"]==field_scalar_data(P,scalar),"word certificate central lift scalar is corrupt")
    else
        ensure(!isnothing(scalar_sign(scalar)),"legacy word certificate cannot encode this central scalar")
    end
    ensure(data["SL_lift_sign"]==scalar_sign(scalar),"word certificate central lift sign is corrupt")
    true
end

"""A verified (not canonical or shortest) GL/SL/PSL word and optional path certificate."""
function matrix_to_word(R::Union{ArithmeticGroupResult,GLGroupResult},g;
                        target=R isa GLGroupResult ? :GL : :SL,simplified=true,certificate=true,
                        cache_path=nothing,max_steps=typemax(Int),timeout=Inf)
    max_steps isa Integer && max_steps>=0 && timeout isa Real && timeout>=0 ||
        throw(WordReductionError(:invalid_input,"nonnegative search budgets are required",nothing))
    g=word_input(R,g,target);P=R.enumeration.problem;alphabet=word_alphabet(R,target,simplified)
    if !isnothing(cache_path) && isfile(cache_path)
        data=JSON3.read(read(cache_path,String),Dict{String,Any})
        try verify_word_certificate(R,g,data) catch err
            throw(WordReductionError(:invalid_certificate,sprint(showerror,err),String(cache_path)))
        end
        data["cache_replay"]=true
        return certificate ? data : Int.(data["word"])
    end
    checkpoint=isnothing(cache_path) ? nothing : String(cache_path)*".partial"
    path=try
        matrix_to_word(R.gl,g;max_cones=max_steps,timeout,checkpoint)
    catch err
        message=sprint(showerror,err)
        kind=occursin("explicit",message) && occursin("limit",message) ? :resource_exhausted : :reduction_failure
        throw(WordReductionError(kind,message,checkpoint))
    end
    word=translated_word(R,path.word,target,simplified)
    value=evaluate_word(alphabet.matrices,word)
    scalar=target==:PSL ? projective_scalar(value,g) : (value==g ? one(P.field) : nothing)
    ensure(!isnothing(scalar),"translated word does not evaluate to the input matrix modulo the requested centre")
    data=JSON3.read(JSON3.write((;schema="av_constructive_word_v2",status="complete",source_hash=R.context.source_hash,
        target=string(target),simplified,alphabet_fingerprint=alphabet.fingerprint,input_matrix=field_matrix_data(P,g),
        word,SL_lift_scalar=field_scalar_data(P,scalar),SL_lift_sign=scalar_sign(scalar),raw_GL_path=path,cache_replay=false,
        guarantee="exact input matrix for GL/SL; input SL lift modulo the full scalar centre for PSL; no canonical or shortest word guarantee")),Dict{String,Any})
    verify_word_certificate(R,g,data)
    isnothing(cache_path) || atomic_json(cache_path,data)
    certificate ? data : word
end
