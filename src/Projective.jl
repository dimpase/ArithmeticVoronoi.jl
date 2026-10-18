export psl_presentation, psl_presentation_data

function psl_presentation(S::SLPresentationCertificate)
    P=S.gl.context.enumeration.problem;I=identity_matrix(P.field,P.rank)
    ensure(evaluate_word(S.generator_matrices,S.minus_identity_word)==-I,"cannot form PSL without the exact central -I word")
    relations=[S.relators;[S.minus_identity_word]]
    signs=[evaluate_word(S.generator_matrices,w)==I ? 1 : -1 for w in relations]
    ensure(all(evaluate_word(S.generator_matrices,w)==sign*I for (w,sign) in zip(relations,signs)),"PSL relator is not scalar +/-I")
    (;generator_matrices=S.generator_matrices,relators=relations,relator_signs=signs,
      killed_minus_identity_word=S.minus_identity_word,source_hash=S.gl.context.source_hash)
end

function psl_presentation_data(P,Q)
    (;schema=PRESENTATION_SCHEMA,status="certified_PSL",group="PSL(2,O) = SL(2,O)/{I,-I}",
      source_hash=Q.source_hash,projective=true,matrix_semantics="determinant-one SL lifts; relations hold modulo +/-I",
      generators=[(;id=i,matrix=field_matrix_data(P,g)) for (i,g) in enumerate(Q.generator_matrices)],
      relators=[(;word=w,SL_lift_sign=Q.relator_signs[i]) for (i,w) in enumerate(Q.relators)],
      quotient_from="certified_SL",killed_minus_identity_word=Q.killed_minus_identity_word,
      proof="quotient of the faithful SL presentation by its identified central order-two element")
end
