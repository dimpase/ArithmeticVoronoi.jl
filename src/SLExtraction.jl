export SLPresentationCertificate, extract_sl_presentation, sl_presentation_data,
       verify_sl_presentation

function verified_units(P)
    ensure(P.degree==2,"unit verification requires a quadratic field")
    U,map=unit_group(P.order;GRH=false)
    if is_hermitian(P)
        m=P.discriminant == -3 ? 6 : (P.discriminant == -4 ? 4 : 2)
        z=m>2 ? P.generator : -one(P.field)
        ensure(order(U)==m && z^m==1 && length(unique(z^k for k in 0:m-1))==m,"incorrect imaginary-quadratic unit group")
        return (;group=U,map,rank=0,torsion_order=m,torsion_generator=z,
            epsilon=one(P.field),epsilon_exponent=BigInt(0),epsilon_torsion=0)
    end
    ensure(ngens(U)==2 && P.field(map(U[1]))==-1,"unexpected real-quadratic unit-group basis")
    generator=P.field(map(U[2]))
    # Among +/-u and +/-u^-1, positive trace and positive w coefficient select
    # the unique unit >1 at w=(root_sum+sqrt(D))/2. This exact normalization also
    # preserves the historical unit choices for D5, D8, D12, D13 and D40.
    candidates=(generator,-generator,inv(generator),-inv(generator))
    epsilon=only(u for u in candidates if tr(u)>0 && coeff(u,1)>0)
    ensure(epsilon in P.order && abs(norm(epsilon))==1,"fundamental unit is not integral with unit norm")
    coordinates=preimage(map,P.order(epsilon))
    exponent=BigInt(coordinates[2]);torsion=Int(mod(coordinates[1],2))
    ensure(abs(exponent)==1,"chosen epsilon is not a generator of the full infinite unit part")
    (;group=U,map,rank=1,torsion_order=2,torsion_generator=-one(P.field),epsilon,epsilon_exponent=exponent,epsilon_torsion=torsion)
end

function finite_unit_exponent(units,delta)
    k=findfirst(k->units.torsion_generator^k==delta,0:units.torsion_order-1)
    ensure(!isnothing(k),"determinant is outside the verified finite unit group")
    k-1
end

function unit_coordinates(P,units,delta)
    coefficients=preimage(units.map,P.order(delta))
    exponent=BigInt(coefficients[2])*units.epsilon_exponent
    sign=Int(mod(coefficients[1]-units.epsilon_torsion*exponent,2))
    ensure(P.field((-1)^sign)*units.epsilon^exponent==delta,"exact determinant unit decomposition failed")
    (;sign,exponent,mod_squares=sign+2Int(mod(exponent,2)))
end

struct SLPresentationCertificate
    gl::PresentationCertificate
    units::Any
    scalar_word::Vector{Int}
    scalar_reduction::Any
    determinant_images::Vector{Int}
    transversal_words::Vector{Vector{Int}}
    schreier_ids::Matrix{Int}
    generator_lifts::Vector{Vector{Int}}
    scalar_corrections::Vector{BigInt}
    generator_matrices::Vector{Any}
    relators::Vector{Vector{Int}}
    relator_provenance::Vector{Any}
    minus_identity_word::Vector{Int}
end

determinant_coset(q,image,sign,modulus)=modulus==0 ? xor(q,image) : mod(q+sign*image,modulus)

function schreier_rewrite(word,images,ids;initial_coset=0,modulus=0)
    result=Int[];coset=initial_coset
    for letter in word
        k=abs(letter);next=determinant_coset(coset,images[k],sign(letter),modulus)
        id=letter>0 ? ids[coset+1,k] : ids[next+1,k]
        id==0 || push!(result,letter>0 ? id : -id)
        coset=next
    end
    (;word=freely_reduce(result),coset)
end
schreier_rewrite(S::SLPresentationCertificate,word)=schreier_rewrite(word,S.determinant_images,S.schreier_ids;
    modulus=S.units.rank==0 ? S.units.torsion_order : 0)

function extract_sl_presentation(C::PresentationCertificate;max_cones=10000,timeout=120.0)
    verify_presentation(C)
    P=C.context.enumeration.problem;units=verified_units(P);epsilon=units.epsilon
    finite=units.rank==0;index=finite ? units.torsion_order : 4
    modulus=finite ? index : 0
    scalar_reduction=finite ? nothing : matrix_to_word(C,epsilon*identity_matrix(P.field,P.rank);max_cones,timeout)
    scalar_word=finite ? Int[] : scalar_reduction.word
    ensure(evaluate_word(C.generator_matrices,scalar_word)==epsilon*identity_matrix(P.field,P.rank),"fundamental scalar word is incorrect")
    images=[finite ? finite_unit_exponent(units,det(g)) : unit_coordinates(P,units,det(g)).mod_squares for g in C.generator_matrices]
    # Prefix-closed Schreier transversal for the actual determinant quotient.
    transversal=Union{Nothing,Vector{Int}}[nothing for _ in 1:index];transversal[1]=Int[]
    queue=[0];cursor=1
    while cursor<=length(queue)
        q=queue[cursor]
        for (j,image) in enumerate(images)
            next=determinant_coset(q,image,1,modulus)
            if isnothing(transversal[next+1])
                transversal[next+1]=[transversal[q+1];j];push!(queue,next)
            end
        end
        cursor+=1
    end
    ensure(length(queue)==index,"determinant images do not generate the full required quotient")
    words=Vector{Int}[x for x in transversal]
    # For imaginary fields use the full finite determinant kernel directly;
    # quotienting by roots of unity here would discard part of SL's centre.
    quotient_relators=finite ? C.relators : [C.relators;[scalar_word]]
    for word in quotient_relators
        ensure(foldl((q,k)->determinant_coset(q,images[abs(k)],sign(k),modulus),word;init=0)==0,"GL/scalar relator is incompatible with the determinant quotient")
    end
    ids=zeros(Int,index,length(images));lifts=Vector{Int}[];corrections=BigInt[];matrices=Any[]
    for q in 0:index-1,j in eachindex(images)
        target=determinant_coset(q,images[j],1,modulus);lift=freely_reduce([words[q+1];j;inverse_word(words[target+1])])
        isempty(lift) && continue
        old=evaluate_word(C.generator_matrices,lift)
        correction=BigInt(0)
        if !finite
            unit=unit_coordinates(P,units,det(old))
            ensure(unit.sign==0 && iseven(unit.exponent),"Schreier generator lift has nonsquare determinant")
            correction=-div(unit.exponent,2)
        end
        matrix=epsilon^correction*old
        ensure(det(matrix)==1,"scalar correction did not produce determinant one")
        check_group_matrix(P,encode_matrix(P,matrix))
        push!(lifts,lift);push!(corrections,correction);push!(matrices,matrix);ids[q+1,j]=length(lifts)
    end
    relators=Vector{Int}[];provenance=Any[]
    for (r,word) in enumerate(quotient_relators),q in 0:index-1
        conjugated=[words[q+1];word;inverse_word(words[q+1])]
        rewrite=schreier_rewrite(conjugated,images,ids;modulus)
        ensure(rewrite.coset==0,"Schreier relation does not return to the identity coset")
        ensure(evaluate_word(matrices,rewrite.word)==identity_matrix(P.field,P.rank),"rewritten SL relation fails exact corrected matrix evaluation")
        isempty(rewrite.word) && continue
        push!(relators,rewrite.word)
        push!(provenance,(;kind=r<=length(C.relators) ? "rewritten_GL_relator" : "rewritten_central_scalar_relator",
            GL_relator=r,coset=q,conjugated_GL_word=conjugated))
    end
    minus_identity=-identity_matrix(P.field,P.rank)
    minus_word=finite_element_word(C,1,minus_identity)
    minus_rewrite=schreier_rewrite(minus_word,images,ids;modulus)
    ensure(minus_rewrite.coset==0 && evaluate_word(matrices,minus_rewrite.word)==minus_identity,"central -I was lost during SL extraction")
    S=SLPresentationCertificate(C,units,scalar_word,scalar_reduction,images,words,ids,lifts,corrections,
        matrices,relators,provenance,minus_rewrite.word)
    verify_sl_presentation(S)
    S
end

function verify_sl_presentation(S::SLPresentationCertificate)
    C=S.gl;P=C.context.enumeration.problem;I=identity_matrix(P.field,P.rank)
    ensure(length(S.transversal_words)==(S.units.rank==0 ? S.units.torsion_order : 4),"wrong determinant quotient index")
    ensure(evaluate_word(C.generator_matrices,S.scalar_word)==S.units.epsilon*I,"corrupt central scalar word")
    for (i,matrix) in enumerate(S.generator_matrices)
        lift=evaluate_word(C.generator_matrices,S.generator_lifts[i])
        ensure(matrix==S.units.epsilon^S.scalar_corrections[i]*lift && det(matrix)==1,"SL generator matrix/lift/correction mismatch")
        check_group_matrix(P,encode_matrix(P,matrix))
    end
    for word in S.relators
        ensure(evaluate_word(S.generator_matrices,word)==I,"SL relator does not evaluate to identity")
    end
    ensure(evaluate_word(S.generator_matrices,S.minus_identity_word)==-I && -I!=I,"SL certificate lost central -I")
    true
end

function sl_presentation_data(S::SLPresentationCertificate)
    P=S.gl.context.enumeration.problem;finite=S.units.rank==0
    (;schema=PRESENTATION_SCHEMA,status="certified_SL",group="SL(2,O), K-linear on the free lattice O^2; central -I retained",
      source_hash=S.gl.context.source_hash,action=P.action_convention,
      method=finite ? "exact Reidemeister-Schreier kernel of GL(2,O) -> O^*; all finite scalar units retained" : "central quotient by the fundamental scalar, then exact index-four Reidemeister-Schreier rewriting and determinant-one rescaling",
      fundamental_unit_coefficients=finite ? nothing : [fraction_string(coeff(S.units.epsilon,k)) for k in 0:1],
      unit_group_proof=finite ? "Hecke unit_group(O;GRH=false); exact order and enumeration of all powers of the torsion generator" : "Hecke unit_group(O;GRH=false); chosen epsilon has free exponent +/-1",
      torsion_unit_order=S.units.torsion_order,
      torsion_unit_generator=[fraction_string(coeff(S.units.torsion_generator,k)) for k in 0:1],
      scalar_GL_word=S.scalar_word,scalar_word_search=S.scalar_reduction,
      GL_generator_determinant_mod_squares=finite ? nothing : S.determinant_images,
      GL_generator_determinant_exponents=finite ? S.determinant_images : nothing,
      determinant_quotient_encoding=finite ? "torsion generator exponent modulo its order; group operation addition" : "torsion sign bit + 2*epsilon exponent parity bit; group operation xor",
      quotient_index=length(S.transversal_words),transversal_GL_words=S.transversal_words,
      schreier_ids=[collect(row) for row in eachrow(S.schreier_ids)],
      generators=[(;id=i,GL_lift_word=S.generator_lifts[i],epsilon_scalar_exponent=string(S.scalar_corrections[i]),
          matrix=field_matrix_data(P,matrix)) for (i,matrix) in enumerate(S.generator_matrices)],
      relators=[(;word,provenance=S.relator_provenance[i]) for (i,word) in enumerate(S.relators)],
      minus_identity_word=S.minus_identity_word,all_determinants_one=true,all_relators_verified=true,
      simplification="raw Schreier presentation retained; no generator elimination yet")
end
