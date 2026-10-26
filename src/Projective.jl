export psl_presentation, psl_presentation_data

# Return the exact scalar relating two determinant-one lifts, or nothing.
# The rank condition excludes scalars which are central in GL but not in SL.
function projective_scalar(A,B)
    size(A)==size(B) || return nothing
    for j in 1:ncols(B),i in 1:nrows(B)
        iszero(B[i,j]) && continue
        scalar=A[i,j]/B[i,j]
        return scalar^nrows(B)==1 && A==scalar*B ? scalar : nothing
    end
    nothing
end
scalar_sign(s)=s==1 ? 1 : s==-1 ? -1 : nothing

function psl_presentation(S::SLPresentationCertificate)
    verify_sl_presentation(S)
    P=S.gl.context.enumeration.problem;I=identity_matrix(P.field,P.rank)
    relations=S.center_order==1 ? copy(S.relators) : [S.relators;[S.center_word]]
    scalars=[projective_scalar(evaluate_word(S.generator_matrices,w),I) for w in relations]
    ensure(all(!isnothing,scalars),"PSL relator is not an SL-central scalar")
    (;generator_matrices=S.generator_matrices,relators=relations,relator_scalars=scalars,
      relator_signs=scalar_sign.(scalars),center_order=S.center_order,center_generator=S.center_generator,
      killed_center_word=S.center_word,killed_minus_identity_word=S.minus_identity_word,
      source_hash=S.gl.context.source_hash)
end

function psl_presentation_data(P,Q)
    (;schema=PRESENTATION_SCHEMA,status="certified_PSL",group="PSL($(P.rank),O) = SL($(P.rank),O)/{zeta*I : zeta^$(P.rank)=1}",
      rank=P.rank,source_hash=Q.source_hash,projective=true,
      matrix_semantics="determinant-one SL lifts; relations hold modulo the full scalar centre",
      generators=[(;id=i,matrix=field_matrix_data(P,g)) for (i,g) in enumerate(Q.generator_matrices)],
      relators=[(;word=w,SL_lift_scalar=field_scalar_data(P,Q.relator_scalars[i]),SL_lift_sign=Q.relator_signs[i]) for (i,w) in enumerate(Q.relators)],
      quotient_from="certified_SL",center_order=Q.center_order,center_generator=field_scalar_data(P,Q.center_generator),
      killed_center_word=Q.killed_center_word,killed_minus_identity_word=Q.killed_minus_identity_word,
      proof="quotient of the faithful SL presentation by its full verified scalar centre")
end
