export field_metadata, abelianization_invariants

function field_metadata(P)
    units=verified_units(P)
    if P.degree==1
        return (;discriminant=1,polynomial_coefficients=[0,1],integral_basis=["1"],
            order_is_maximal=true,basis_discriminant=1,fundamental_unit=nothing,
            fundamental_unit_norm=nothing,unit_rank=0,fundamental_unit_free_exponent=nothing,
            fundamental_unit_torsion=nothing,torsion_unit_order=2,torsion_unit_generator=["-1"],
            ordinary_class_number="1",ordinary_class_group_invariants=String[],
            class_group_kind="ordinary",narrow_class_group=nothing,
            proof_settings="Z is a principal ideal domain with units +/-1",
            lmfdb_label=nothing,lmfdb_generator_map=nothing,lmfdb_status="not_requested")
    end
    C,map=class_group(P.order;GRH=false)
    h=BigInt(order(C))
    ensure(h>0,"ordinary class-group computation did not return a finite group")
    (;discriminant=P.discriminant,polynomial_coefficients=[P.root_product,-P.root_sum,1],
      integral_basis=["1","w"],order_is_maximal=discriminant(P.order)==P.discriminant,
      basis_discriminant=P.root_sum^2-4P.root_product,
      fundamental_unit=units.rank==0 ? nothing : [fraction_string(coeff(units.epsilon,k)) for k in 0:1],
      fundamental_unit_norm=units.rank==0 ? nothing : fraction_string(norm(units.epsilon)),unit_rank=units.rank,
      fundamental_unit_free_exponent=units.rank==0 ? nothing : string(units.epsilon_exponent),
      fundamental_unit_torsion=units.rank==0 ? nothing : units.epsilon_torsion,
      torsion_unit_order=units.torsion_order,
      torsion_unit_generator=[fraction_string(coeff(units.torsion_generator,k)) for k in 0:1],
      ordinary_class_number=string(h),
      ordinary_class_group_invariants=[string(order(C[i])) for i in 1:ngens(C) if order(C[i])>1],
      class_group_kind="ordinary",narrow_class_group=nothing,
      proof_settings="Hecke class_group(O; GRH=false), unit_group(O; GRH=false); exact maximal-order discriminant",
      lmfdb_label=nothing,lmfdb_generator_map=nothing,
      lmfdb_status="not_requested")
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
