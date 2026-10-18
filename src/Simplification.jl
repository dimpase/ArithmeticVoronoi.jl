export simplify_matrix_presentation, verify_simplification

"""Tietze simplification with explicit old/new generator words and matrix maps."""
function simplify_matrix_presentation(matrices,relations;label="GL",projective=false)
    F=free_group(length(matrices))
    G,projection=quo(F,[F(word) for word in relations])
    H,isomorphism=simplified_fp_group(G)
    ensure(is_bijective(isomorphism),"Tietze map is not a certified group isomorphism")
    forward=[letters(isomorphism(g)) for g in gens(G)]
    backward=[letters(preimage(isomorphism,h)) for h in gens(H)]
    images=Any[evaluate_word(matrices,word) for word in backward]
    relators_new=[letters(word) for word in relators(H)]
    certificate=(;label,projective,algorithm="OSCAR simplified_fp_group / GAP IsomorphismSimplifiedFpGroup; Tietze isomorphism",
        isomorphism_verified=true,old_generator_count=length(matrices),old_relator_count=length(relations),
        forward_generator_words=forward,backward_generator_words=backward,
        generator_matrices=images,relators=relators_new,
        gap_version=string(Oscar.GAP.gap_to_julia(Oscar.GAP.Globals.GAPInfo.Version)))
    verify_simplification(matrices,relations,certificate)
    certificate
end

function substitute_word(word,images)
    result=Int[]
    for letter in word
        append!(result,letter>0 ? images[letter] : inverse_word(images[-letter]))
    end
    freely_reduce(result)
end

function verify_simplification(old_matrices,old_relations,simplified)
    ensure(simplified.isomorphism_verified,"simplification lacks its Tietze isomorphism")
    new=simplified.generator_matrices;I=identity_matrix(base_ring(first(old_matrices)),nrows(first(old_matrices)))
    projective=hasproperty(simplified,:projective) && simplified.projective
    equal_image(A,B)=A==B || (projective && A==-B)
    ensure(length(simplified.forward_generator_words)==length(old_matrices),"incomplete old-to-new generator map")
    ensure(length(simplified.backward_generator_words)==length(new),"incomplete new-to-old generator map")
    for (i,word) in enumerate(simplified.forward_generator_words)
        ensure(equal_image(evaluate_word(new,word),old_matrices[i]),"old-to-new matrix map does not commute")
    end
    for (i,word) in enumerate(simplified.backward_generator_words)
        ensure(equal_image(evaluate_word(old_matrices,word),new[i]),"new-to-old matrix map does not commute")
    end
    for word in simplified.relators;ensure(equal_image(evaluate_word(new,word),I),"simplified relator fails exact matrix evaluation");end
    for word in old_relations
        ensure(equal_image(evaluate_word(new,substitute_word(word,simplified.forward_generator_words)),I),"rewritten old relation fails the simplified matrix map")
    end
    true
end

function simplification_data(P,S)
    (;schema=PRESENTATION_SCHEMA,status="certified_Tietze_isomorphism",group=S.label,algorithm=S.algorithm,
      projective=hasproperty(S,:projective) && S.projective,
      old_generator_count=S.old_generator_count,old_relator_count=S.old_relator_count,
      generators=[(;id=i,matrix=field_matrix_data(P,g),word_in_old_generators=S.backward_generator_words[i]) for (i,g) in enumerate(S.generator_matrices)],
      old_generators_as_words_in_new=S.forward_generator_words,relators=S.relators,
      isomorphism_verified=S.isomorphism_verified,matrix_maps_verified=true,gap_version=S.gap_version)
end
