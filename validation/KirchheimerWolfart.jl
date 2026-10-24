# Optional reference comparison tooling. Never included by ArithmeticVoronoi.
module KirchheimerWolfart
using ArithmeticVoronoi, Oscar, JSON3
import ArithmeticVoronoi: ensure, evaluate_word, inverse_word, freely_reduce,
    substitute_word, verified_units, unit_coordinates, fraction_string,
    field_matrix_data, atomic_json, PACKAGE_ROOT, file_sha, check_projective_comparison

export kw_reference, compare_kw, reference_matrix_word

function letter_power(letter,n;limit=1000000)
    abs(n)<=limit || throw(ArgumentError("expanded word exceeds explicit length limit"))
    fill(n<0 ? -letter : letter,Int(abs(n)))
end

function kw_reference(P)
    coefficients=Dict(5=>(1,1),8=>(1,2),12=>(-1,4),13=>(1,3))
    ensure(haskey(coefficients,P.discriminant),"Kirchheimer-Wolfart reference applies only to D5,D8,D12,D13")
    alpha,beta=coefficients[P.discriminant];e=verified_units(P).epsilon;K=P.field;I=identity_matrix(K,2)
    ensure(e^2==alpha+beta*e && e^3==alpha*beta+(alpha+beta^2)*e,"reference unit coefficient identity failed")
    U=matrix(K,[1 1;0 1]);V=matrix(K,[1 e;0 1]);H=diagonal_matrix(K,[e,inv(e)]);J=matrix(K,[0 1;-1 0])
    matrices=Any[U,V,H,J,-I]
    common=[[1,2,-1,-2],vcat([3,1,-3],letter_power(1,-alpha),letter_power(2,-beta)),
        vcat([3,2,-3],letter_power(1,-alpha*beta),letter_power(2,-(alpha+beta^2)))]
    psl=[common;[[4,4],[3,4,3,4],repeat([4,1],3),repeat([3,4,2],3)]]
    sl=[common;[[5,5]];[[5,j,-5,-j] for j in 1:4];[[4,4,-5],[3,4,3,4,-5],repeat([4,1],3),repeat([3,4,2],3)]]
    if P.discriminant==12
        extra=repeat([4,2,-1,4,2,-1,-1,-1],2)
        push!(psl,extra);push!(sl,[extra;-5])
    end
    ensure(all(evaluate_word(matrices,w)==I for w in sl),"reference SL relator has an incorrect central sign")
    ensure(all(evaluate_word(matrices,w) in (I,-I) for w in psl),"reference PSL relator is not projectively trivial")
    (;generator_matrices=matrices,SL_relators=sl,PSL_relators=psl,alpha,beta,
      generator_names=["U=tau(1)","V=tau(e)","H=diag(e,e^-1)","J=[0,1;-1,0]","Z=-I"],
      epsilon=e,proof="published PSL theorem (9.1)/(10.1)/(11.3)/(12.2), with the exact central extension described in section 13")
end

"""A finite exact row reduction for the four reference fields; never used by the Voronoi engine."""
function reference_matrix_word(P,g;max_steps=1000,division_radius=6)
    ensure(P.discriminant in (5,8,12,13),"reference row reduction is restricted to the four comparison fields")
    check_group_matrix(P,encode_matrix(P,g));ensure(det(g)==1,"reference input must have determinant one")
    ref=kw_reference(P);units=verified_units(P);K=P.field;I=identity_matrix(K,2)
    function translation(t)
        a,b=coeff(t,0),coeff(t,1);c=coeff(units.epsilon,0)
        ensure(denominator(a)==1 && denominator(b)==1,"nonintegral reference translation")
        [letter_power(1,BigInt(ZZ(a-b*c)));letter_power(2,BigInt(ZZ(b)))]
    end
    M=deepcopy(g);left=Int[];steps=Any[]
    for number in 0:max_steps
        if M[2,1]==0
            unit=unit_coordinates(P,units,M[1,1])
            terminal=[letter_power(5,unit.sign);letter_power(3,unit.exponent);translation(M[1,2]/M[1,1])]
            result=freely_reduce([inverse_word(left);terminal])
            ensure(evaluate_word(ref.generator_matrices,result)==g,"reference row-reduction word has the wrong matrix")
            return (;word=result,steps,method="exact norm-decreasing row operations; finite input-specific division witnesses")
        end
        number<max_steps || throw(ArithmeticValidationError("reference row reduction reached its explicit step limit"))
        a,b=M[1,1],M[2,1];x=a/b
        floors=[fld(BigInt(numerator(coeff(x,k))),BigInt(denominator(coeff(x,k)))) for k in 0:1]
        q=nothing;remainder=nothing
        for radius in 0:division_radius
            for j in -radius:radius+1,i in -radius:radius+1
                trial=K(floors[1]+i)+(floors[2]+j)*P.generator;r=a-trial*b
                if abs(norm(r))<abs(norm(b));q=trial;remainder=r;break;end
            end
            isnothing(q) || break
        end
        isnothing(q) && throw(ArithmeticValidationError("reference division found no decreasing quotient within its explicit search radius"))
        factor=[4;translation(-q)]
        left=freely_reduce([factor;left]);M=ref.generator_matrices[4]*matrix(K,[1 -q;0 1])*M
        ensure(M[2,1]==-remainder,"reference Euclidean row convention failed")
        push!(steps,(;quotient=[fraction_string(coeff(q,k)) for k in 0:1],
            previous_lower_norm=fraction_string(abs(norm(b))),next_lower_norm=fraction_string(abs(norm(remainder)))))
    end
    error("unreachable")
end

function kw_source_hashes(source_directory=nothing)
    record=JSON3.read(read(joinpath(PACKAGE_ROOT,"references/kw_sources.json"),String),Dict{String,Any})
    hashes=record["sha256"]
    if !isnothing(source_directory)
        for (name,expected) in hashes
            path=joinpath(source_directory,name)
            ensure(isfile(path),"missing reference transcription: $path")
            ensure(file_sha(path)==expected,"reference transcription hash differs: $name")
        end
    end
    hashes
end

function compare_kw(R;directory=joinpath(R.directory,"KW"),reference_source_directory=nothing)
    source_hashes=kw_source_hashes(reference_source_directory)
    mkpath(directory);P=R.enumeration.problem;ref=kw_reference(P);small=R.simplified[:SL]
    forward=Any[]
    for (i,g) in enumerate(ref.generator_matrices)
        path=joinpath(directory,"reference_generator_$i.json")
        push!(forward,matrix_to_word(R,g;target=:SL,cache_path=path))
    end
    backward=[reference_matrix_word(P,g) for g in small.generator_matrices]
    A=[Int.(x["word"]) for x in forward];B=[x.word for x in backward]
    I=identity_matrix(P.field,2)
    for (i,g) in enumerate(ref.generator_matrices)
        ensure(evaluate_word(small.generator_matrices,A[i])==g,"reference-to-computed matrix map failed")
        ensure(evaluate_word(ref.generator_matrices,substitute_word(A[i],B))==g,"reference round trip changed its matrix")
    end
    for (i,g) in enumerate(small.generator_matrices)
        ensure(evaluate_word(ref.generator_matrices,B[i])==g,"computed-to-reference matrix map failed")
        ensure(evaluate_word(small.generator_matrices,substitute_word(B[i],A))==g,"computed round trip changed its matrix")
    end
    for w in ref.SL_relators
        ensure(evaluate_word(small.generator_matrices,substitute_word(w,A))==I,"substituted reference relator failed")
    end
    for w in small.relators
        ensure(evaluate_word(ref.generator_matrices,substitute_word(w,B))==I,"substituted computed relator failed")
    end
    abref=abelianization_invariants(5,ref.SL_relators);abcomputed=abelianization_invariants(length(small.generator_matrices),small.relators)
    ensure(abref.free_rank==abcomputed.free_rank && abref.torsion==abcomputed.torsion,"reference/computed SL abelianizations differ")
    pslwords=[matrix_to_word(R,g;target=:PSL)["word"] for g in ref.generator_matrices[1:4]]
    pslbackward=[reference_matrix_word(P,g) for g in R.simplified[:PSL].generator_matrices]
    pslA=[Int.(w) for w in pslwords]
    pslB=[freely_reduce([k for k in b.word if abs(k)!=5]) for b in pslbackward]
    check_projective_comparison(ref.generator_matrices[1:4],R.simplified[:PSL].generator_matrices,
        pslA,pslB,ref.PSL_relators,R.simplified[:PSL].relators)
    data=(;schema="av_kw_comparison_v1",discriminant=P.discriminant,status="certified_by_independent_faithfulness",
        source_hash=R.context.source_hash,reference_generator_names=ref.generator_names,
        reference_matrices=field_matrix_data.(Ref(P),ref.generator_matrices),SL_reference_relators=ref.SL_relators,
        PSL_reference_relators=ref.PSL_relators,reference_to_computed_SL_words=A,computed_SL_to_reference_words=B,
        computed_SL_to_reference_reductions=backward,reference_to_computed_PSL_words=pslwords,
        computed_PSL_to_reference_words=pslB,
        PSL_reference_generator_names=ref.generator_names[1:4],
        computed_PSL_to_reference_SL_lift_words=[b.word for b in pslbackward],
        computed_PSL_to_reference_lift_signs=[evaluate_word(ref.generator_matrices[1:4],w)==g ? 1 : -1 for (w,g) in zip(pslB,R.simplified[:PSL].generator_matrices)],
        field_isomorphism="identity on the stated w basis; reference e has the coefficients recorded in field.json",
        conjugation="none",J_sign="J is [0,1;-1,0], the negative of [0,-1;1,0]",
        SL_abelianization=abcomputed,
        proof_method="both presentations independently faithful in the same SL matrices: geometric coverage plus group-action theorem for the computation; published PSL theorem and central-kernel argument for the reference. Exact commuting matrix substitutions therefore give mutually inverse group maps. No formal relator derivations or Tietze chain between the two reference alphabets is claimed.",
        reference_source_hashes=source_hashes,
        reference_source_verification=isnothing(reference_source_directory) ? "recorded_hashes_only" : "checked_local_sources",
        D8_transcription_note="English (10.1) has final exponent -2; German transcription has -5. Generate -5 from e^3=2+5e. The erroneous word is tau(3e). Original scan not yet inspected; no published-paper error is asserted.")
    atomic_json(joinpath(directory,"comparison.json"),data)
    atomic_json(joinpath(directory,"reference.json"),(;schema="av_kw_reference_v1",discriminant=P.discriminant,
        names=ref.generator_names,matrices=field_matrix_data.(Ref(P),ref.generator_matrices),SL_relators=ref.SL_relators,
        PSL_relators=ref.PSL_relators,alpha=ref.alpha,beta=ref.beta,proof=ref.proof))
    R.run["status"]["literature_comparison"]="certified_by_independent_faithfulness"
    atomic_json(joinpath(R.directory,"run.json"),R.run)
    println("D$(P.discriminant): Kirchheimer-Wolfart maps verified in both directions.");flush(stdout)
    data
end

end
