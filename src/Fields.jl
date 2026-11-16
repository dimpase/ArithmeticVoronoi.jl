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

# An exponent-sum row with a coefficient +/-1 eliminates that generator by
# integral row operations. Sparse elimination removes redundant presentation
# generators before the dense Smith form, without changing the presented group.
function abelianization_invariants(generator_count,relations)
    generator_count isa Integer && generator_count>=0 ||
        throw(ArgumentError("generator_count must be nonnegative"))
    rows=Dict{Int,BigInt}[];seen=Set{Any}()
    for word in relations
        row=Dict{Int,BigInt}()
        for letter in word
            j=abs(letter)
            1<=j<=generator_count || throw(ArgumentError("relator letter is outside the generator range"))
            row[j]=get(row,j,big(0))+sign(letter)
        end
        filter!(p->!iszero(last(p)),row)
        isempty(row) && continue
        key=Tuple(sort!(collect(row);by=first))
        orientation=sign(last(first(key)))
        normalized=Tuple((j,orientation*a) for (j,a) in key)
        normalized in seen && continue
        push!(seen,normalized);push!(rows,row)
    end
    incoming=[Set{Int}() for _ in 1:generator_count]
    for (i,row) in enumerate(rows),j in keys(row);push!(incoming[j],i);end
    active=trues(generator_count)
    sizes=zeros(Int,length(rows));buckets=Dict{Int,Set{Int}}();minimum_size=Ref(typemax(Int))
    function enqueue(i)
        old=sizes[i];old>0 && delete!(buckets[old],i)
        size=any(x->abs(x)==1,values(rows[i])) ? length(rows[i]) : 0
        sizes[i]=size
        if size>0
            push!(get!(Set{Int},buckets,size),i);minimum_size[]=min(minimum_size[],size)
        end
    end
    foreach(enqueue,eachindex(rows))
    while any(!isempty,values(buckets))
        while !haskey(buckets,minimum_size[]) || isempty(buckets[minimum_size[]]);minimum_size[]+=1;end
        i=first(buckets[minimum_size[]]);delete!(buckets[minimum_size[]],i);sizes[i]=0
        pivot=rows[i]
        columns=[j for (j,a) in pivot if abs(a)==1]
        j=argmin(j->(length(incoming[j]),j),columns);unit=pivot[j]
        for column in keys(pivot);delete!(incoming[column],i);end
        rows[i]=Dict{Int,BigInt}();active[j]=false
        for k in collect(incoming[j])
            row=rows[k];factor=row[j]*unit
            for (column,a) in pivot
                before=get(row,column,big(0));after=before-factor*a
                if iszero(after)
                    delete!(row,column);delete!(incoming[column],k)
                else
                    row[column]=after
                    iszero(before) && push!(incoming[column],k)
                end
            end
            enqueue(k)
        end
        @assert isempty(incoming[j])
    end
    columns=findall(active);positions=Dict(j=>i for (i,j) in enumerate(columns))
    remaining=filter(!isempty,rows)
    A=zero_matrix(ZZ,length(remaining),length(columns))
    for (i,row) in enumerate(remaining),(j,a) in row;A[i,positions[j]]=a;end
    D=snf(A)
    diagonals=[abs(D[i,i]) for i in 1:min(nrows(D),ncols(D)) if D[i,i]!=0]
    (;free_rank=length(columns)-length(diagonals),torsion=[string(x) for x in diagonals if x>1],
      method="exact unit elimination followed by Smith normal form of the relator exponent-sum matrix")
end
