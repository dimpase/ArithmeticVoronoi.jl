export field_metadata, abelianization_invariants

function field_metadata(P)
    units=verified_units(P)
    C,map=class_group(P.order;GRH=false)
    h=BigInt(order(C))
    ensure(h==(P.discriminant==40 ? 2 : 1),"unexpected ordinary class number; inspect the exact class-group computation")
    (;discriminant=P.discriminant,polynomial_coefficients=[P.root_product,-P.root_sum,1],
      integral_basis=["1","w"],order_is_maximal=discriminant(P.order)==P.discriminant,
      basis_discriminant=P.root_sum^2-4P.root_product,
      fundamental_unit=[fraction_string(coeff(units.epsilon,k)) for k in 0:1],
      fundamental_unit_norm=fraction_string(norm(units.epsilon)),unit_rank=1,
      fundamental_unit_free_exponent=string(units.epsilon_exponent),
      fundamental_unit_torsion=units.epsilon_torsion,
      ordinary_class_number=string(h),
      ordinary_class_group_invariants=[string(order(C[i])) for i in 1:ngens(C) if order(C[i])>1],
      class_group_kind="ordinary",narrow_class_group=nothing,
      proof_settings="Hecke class_group(O; GRH=false), unit_group(O; GRH=false); exact maximal-order discriminant",
      lmfdb_label=nothing,lmfdb_generator_map=nothing,
      lmfdb_status="unverified; LMFDB returned browser challenges on the D5/D8 verification probes")
end

function abelianization_invariants(generator_count,relations)
    A=zero_matrix(ZZ,length(relations),generator_count)
    for (i,word) in enumerate(relations),letter in word
        A[i,abs(letter)]+=sign(letter)
    end
    D=snf(A)
    diagonals=[abs(D[i,i]) for i in 1:min(nrows(D),ncols(D)) if D[i,i]!=0]
    (;free_rank=generator_count-length(diagonals),torsion=[string(x) for x in diagonals if x>1],
      method="Smith normal form of the relator exponent-sum matrix")
end
