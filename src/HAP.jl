export VoronoiCell, VoronoiComplex, voronoi_complex, validate_voronoi_complex,
       hap_complex, hap_voronoi_complex, enable_hap_julia_backend

import Oscar: GAP

# Bounded faces of the relative Ryshkov polyhedron, dual to the well-rounded
# Voronoi fan. A k-cell is specified by its complete minimum set S, with
# k = dim(form space) - rank(evaluation(S)). Only K-spanning S are retained.
# See Schönnenbeck, "Resolutions for unit groups of orders", arXiv:1609.08835.
struct VoronoiCell
    vectors::Vector{Vector{BigInt}}
    point::Any
    tangent::QQMatrix
    metric::Any
end

struct VoronoiComplex
    enumeration::ArithmeticVoronoiEnumeration
    group::Symbol
    cells::Vector{Vector{VoronoiCell}}
    stabilizers::Vector{Vector{Vector{Any}}}
    # (signed cell index, matrix): a left translate of an oriented boundary cell.
    boundaries::Vector{Vector{Vector{Tuple{Int,Any}}}}
    generators::Vector{Any}
end

cell_dimension(c::VoronoiCell)=ncols(c.tangent)
cell_vector_key(V)=Tuple(sort!(collect(signed_key(v) for v in V)))
cell_image(P,g,V)=[apply_integral(encode_matrix(P,g),v) for v in V]

function cell_basis(P,V)
    B=zero_matrix(P.field,P.rank,0)
    for v in V
        next=hcat(B,matrix(P.field,P.rank,1,field_vector(P,v)))
        rank(next)>ncols(B) || continue
        B=next
        ncols(B)==P.rank && return B
    end
    nothing
end

function make_voronoi_cell(P,V,F)
    vectors=[collect(v) for v in cell_vector_key(V)]
    ensure(!isnothing(cell_basis(P,vectors)),"cell minimum set does not span K^$(P.rank)")
    rays=evaluation_rays(P,vectors)
    _,T=nullspace(rays)
    S=sum(vector_outer(P,field_vector(P,v)) for v in vectors)
    ensure(positive_definite(trace_gram(P,S)),"cell has no positive interior in the dual cone")
    VoronoiCell(vectors,F,T,inv(S))
end

function transform_cell(P,c,g)
    make_voronoi_cell(P,cell_image(P,g,c.vectors),act_form(P,c.point,inv(g)))
end

# An integral linear map is determined by the images of a K-basis. Backtracking
# checks all source vectors as soon as their images are determined. For equal
# minimum sets, their canonical metrics also constrain every pair of images.
# Both signs of every target vector are retained.
function cell_maps(P,source,target; determinant_one=false, first_only=false)
    B=cell_basis(P,source)
    ensure(!isnothing(B),"transporter source is not well rounded")
    inverse=inv(B)
    signed=[field_vector(P,v) for v in [target;[-v for v in target]]]
    W=Set(Tuple(v) for v in signed)
    coordinates=[inverse*matrix(P.field,P.rank,1,field_vector(P,v)) for v in source]
    stages=[[c for c in coordinates if findlast(i->!iszero(c[i,1]),1:P.rank)==j] for j in 1:P.rank]
    Y=zero_matrix(P.field,P.rank,P.rank)
    equal_size=length(source)==length(target)
    inner(u,M,v)=sum(form_conjugate(P,u[i])*M[i,j]*v[j] for i in 1:P.rank,j in 1:P.rank)
    metric=equal_size ? inv(sum(vector_outer(P,field_vector(P,v)) for v in target)) : nothing
    gram=equal_size ? form_adjoint(P,B)*inv(sum(vector_outer(P,field_vector(P,v)) for v in source))*B : nothing
    result=Any[]
    function extend(k)
        for v in signed
            if equal_size
                inner(v,metric,v)==gram[k,k] || continue
                all(inner(Y[:,j],metric,v)==gram[j,k] for j in 1:k-1) || continue
            end
            Y[:,k]=v
            rank(Y[:,1:k])==k || continue
            all(Tuple((Y[:,1:k]*c[1:k,:])[:,1]) in W for c in stages[k]) || continue
            if k<P.rank
                extend(k+1) && return true
            else
                g=Y*inverse
                delta=det(g)
                determinant_one && delta!=1 && continue
                (P.degree==1 ? abs(delta)==1 : abs(norm(delta))==1) || continue
                all(x->isone(denominator(x)),encode_matrix(P,g)) || continue
                all(Tuple((Y*c)[:,1]) in W for c in coordinates) || continue
                push!(result,g)
                first_only && return true
            end
        end
        false
    end
    extend(1)
    result
end

cell_invariant(P,c)=(cell_dimension(c),length(c.vectors),det(trace_gram(P,c.metric)))

function add_cell_orbit!(P,cells,buckets,c;determinant_one=false)
    key=cell_invariant(P,c)
    candidates=get!(buckets,key,VoronoiCell[])
    for r in candidates
        isempty(cell_maps(P,c.vectors,r.vectors;determinant_one,first_only=true)) || return false
    end
    push!(candidates,c)
    push!(cells[cell_dimension(c)+1],c)
    true
end

function native_cell_orbits(E;progress=false)
    P=E.problem
    cells=[VoronoiCell[] for _ in 1:length(P.form_basis)]
    buckets=Dict{Any,Vector{VoronoiCell}}()
    for (i,r) in enumerate(E.representatives)
        edges=[e for e in E.neighbours if e.source==i]
        across=[act_form(P,form_from_trace(P,edges[f.neighbour_orbit].neighbour_gram),
                         f.to_orbit_representative) for f in r.facets]
        queue=[collect(eachindex(r.minimal_vectors))]
        seen=Set{Tuple}([Tuple(only(queue))])
        cursor=1
        while cursor<=length(queue)
            indices=queue[cursor]
            V=r.minimal_vectors[indices]
            if !isnothing(cell_basis(P,V))
                incident=findall(f->issubset(indices,f.indices),r.facets)
                # A strict convex combination of the vertex and all neighbouring
                # vertices across these facets has exactly this minimum set.
                F=isempty(incident) ? r.form :
                    (r.form+sum(across[j] for j in incident)/length(incident))/2
                ensure(findall(v->quadratic_value(trace_gram(P,F),v)==1,r.minimal_vectors)==indices,
                       "cell interior has an incorrect minimum set")
                add_cell_orbit!(P,cells,buckets,make_voronoi_cell(P,V,F))
                for f in r.facets
                    face=intersect(indices,f.indices)
                    isempty(face) && continue
                    key=Tuple(face)
                    key in seen && continue
                    push!(seen,key);push!(queue,face)
                end
            end
            cursor+=1
        end
        if progress
            println("Native Voronoi cells: perfect form ",i,"/",length(E.representatives),
                    "; cell orbits by dimension ",length.(cells));flush(stdout)
        end
    end
    while isempty(last(cells));pop!(cells);end
    cells
end

function cell_tangent_image(P,c,g)
    inverse=inv(g)
    T=c.tangent
    columns=[form_coordinates(P,act_form(P,form_from_trace(P,native_direction(P,T[:,j])),inverse))
             for j in 1:ncols(T)]
    matrix(QQ,nrows(T),ncols(T),[columns[j][i] for i in 1:nrows(T) for j in 1:ncols(T)])
end

function cell_coordinates(c,V)
    T=c.tangent
    X=inv(transpose(T)*T)*transpose(T)*V
    ensure(T*X==V,"cell incidence or action leaves its affine tangent space")
    X
end

function cell_orientation(P,c,g)
    cell_dimension(c)==0 && return 1
    d=det(cell_coordinates(c,cell_tangent_image(P,c,g)))
    ensure(abs(d)==1,"finite cell stabilizer has nonunit orientation determinant")
    Int(sign(d))
end

function incidence_sign(P,upper,lower,g)
    point=act_form(P,lower.point,inv(g))
    outward=form_coordinates(P,point-upper.point)
    tangent=cell_tangent_image(P,lower,g)
    frame=hcat(matrix(QQ,length(outward),1,outward),tangent)
    d=det(cell_coordinates(upper,frame))
    ensure(!iszero(d),"degenerate oriented cell incidence")
    Int(sign(d))
end

"""
    voronoi_complex(E; group=:GL, progress=false)

Construct all bounded Ryshkov cells, their finite stabilizers, and signed
equivariant boundaries from a complete perfect-form enumeration.
Supports rational/quadratic GL and imaginary-quadratic SL on the free lattice O^n.
No Magma, stored complexes, or C++ exporter is used by the native enumeration.
"""
function voronoi_complex(E::ArithmeticVoronoiEnumeration;group=:GL,progress=false)
    P=E.problem
    group in (:GL,:SL) || throw(ArgumentError("choose group=:GL or :SL"))
    group==:SL && !is_hermitian(P) && throw(ArgumentError("SL cell complexes currently require an imaginary quadratic field"))
    ensure(E.completion_status==:complete && isempty(E.pending_facets),"cell complex requires a complete enumeration")
    cells=native_cell_orbits(E;progress)
    generators=Any[g for r in E.representatives for g in r.stabilizer]
    append!(generators,[inv(e.transporter) for e in E.neighbours])
    if group==:SL
        units=verified_units(P)
        reps=[diagonal_matrix(P.field,[units.torsion_generator^j;fill(one(P.field),P.rank-1)]) for j in 0:units.torsion_order-1]
        split=[VoronoiCell[] for _ in cells];buckets=Dict{Any,Vector{VoronoiCell}}()
        for layer in cells,c in layer,t in reps
            add_cell_orbit!(P,split,buckets,transform_cell(P,c,t);determinant_one=true)
        end
        cells=split
        generators=Any[t*g*inv(reps[finite_unit_exponent(units,det(t*g))+1]) for t in reps for g in generators]
        ensure(all(g->det(g)==1,generators),"Schreier generator has nontrivial determinant")
    end
    unique!(generators)
    stabilizers=[[cell_maps(P,c.vectors,c.vectors;determinant_one=group==:SL) for c in layer] for layer in cells]
    boundaries=[[Tuple{Int,Any}[] for _ in layer] for layer in cells]
    for k in 1:length(cells)-1,(i,upper) in enumerate(cells[k+1])
        for (j,lower) in enumerate(cells[k])
            seen=Set{Any}()
            for h in cell_maps(P,upper.vectors,lower.vectors;determinant_one=group==:SL)
                g=inv(h)
                key=cell_vector_key(cell_image(P,g,lower.vectors))
                key in seen && continue
                push!(seen,key)
                push!(boundaries[k+1][i],(incidence_sign(P,upper,lower,g)*j,g))
            end
        end
        ensure(!isempty(boundaries[k+1][i]),"positive-dimensional cell has no boundary")
        if progress
            println("Native Voronoi boundary: dimension ",k,", orbit ",i,"/",length(cells[k+1]));flush(stdout)
        end
    end
    C=VoronoiComplex(E,group,cells,stabilizers,boundaries,generators)
    validate_voronoi_complex(C)
    C
end

"""Check stabilizers, oriented incidences and d²=0 in the integral G-modules."""
function validate_voronoi_complex(C::VoronoiComplex)
    P=C.enumeration.problem
    for g in C.generators
        check_group_matrix(P,encode_matrix(P,g))
        C.group==:SL && ensure(det(g)==1,"SL generator has nontrivial determinant")
    end
    for (k,layer) in enumerate(C.cells),(i,c) in enumerate(layer)
        group=C.stabilizers[k][i]
        identity=identity_matrix(P.field,P.rank)
        ensure(identity in group,"cell stabilizer lost its identity")
        if C.group==:GL || iseven(P.rank)
            ensure(-identity in group,"cell stabilizer lost central -I")
        end
        for g in group
            check_group_matrix(P,encode_matrix(P,g))
            C.group==:SL && ensure(det(g)==1,"SL cell stabilizer has nontrivial determinant")
            ensure(cell_vector_key(cell_image(P,g,c.vectors))==cell_vector_key(c.vectors),"invalid cell stabilizer")
            cell_orientation(P,c,g)
        end
        for (j,g) in C.boundaries[k][i]
            check_group_matrix(P,encode_matrix(P,g))
            C.group==:SL && ensure(det(g)==1,"SL boundary transporter has nontrivial determinant")
            lower=C.cells[k-1][abs(j)]
            ensure(issubset(Set(cell_vector_key(c.vectors)),Set(cell_vector_key(cell_image(P,g,lower.vectors)))),
                   "invalid boundary minimum-set inclusion")
            ensure(sign(j)==incidence_sign(P,c,lower,g),"incorrect boundary orientation")
        end
        k<3 && continue
        coefficients=Dict{Any,Tuple{Any,Int}}()
        for (j,g) in C.boundaries[k][i],(l,h) in C.boundaries[k-1][abs(j)]
            lower=C.cells[k-2][abs(l)];a=g*h
            key=(abs(l),cell_vector_key(cell_image(P,a,lower.vectors)))
            reference,value=get(coefficients,key,(a,0))
            s=inv(reference)*a
            ensure(s in C.stabilizers[k-2][abs(l)],"boundary comparison is not a stabilizer")
            coefficients[key]=(reference,value+sign(j)*sign(l)*cell_orientation(P,lower,s))
        end
        ensure(all(iszero(last(x)) for x in values(coefficients)),"equivariant boundary does not square to zero")
    end
    true
end

function cell_generators(P,elements)
    generators=finite_group_record(P,elements;with_table=false).generators
    # Keep the matrix identity so GAP's trivial group has the correct matrix type.
    isempty(generators) ? Any[identity_matrix(P.field,P.rank)] : generators
end

function hap_complex(C::VoronoiComplex;resolution_length=100)
    resolution_length isa Integer && resolution_length>0 || throw(ArgumentError("resolution_length must be positive"))
    validate_voronoi_complex(C)
    GAP.Packages.load("hap") || error("HAP is not installed; use GAP.Packages.install(\"hap\")")
    GAP.Globals.Read(GAP.GapObj(joinpath(PACKAGE_ROOT,"gap","ArithmeticVoronoi.g")))
    P=C.enumeration.problem
    elts=Any[identity_matrix(P.field,P.rank)]
    indices=Dict(matrix_key(encode_matrix(P,only(elts)))=>1)
    function index(g)
        key=matrix_key(encode_matrix(P,g))
        get!(indices,key) do
            push!(elts,g);length(elts)
        end
    end
    wire(g)=[[BigInt(ZZ(x)) for x in row] for row in matrix_rows(encode_matrix(P,g))]
    stabs=[[cell_generators(P,S) for S in layer] for layer in C.stabilizers]
    evens=[[cell_generators(P,[g for g in C.stabilizers[k][i] if cell_orientation(P,c,g)==1])
            for (i,c) in enumerate(layer)] for (k,layer) in enumerate(C.cells)]
    # HAP may need every stabilizer element in the initial element table.
    for layer in C.stabilizers,S in layer,g in S;index(g);end
    boundaries=[[[[j,index(g)] for (j,g) in b] for b in layer] for layer in C.boundaries]
    data=Dict("dimensions"=>length.(C.cells),"elts"=>wire.(elts),
        "stabilizers"=>[[wire.(S) for S in layer] for layer in stabs],
        "evenStabilizers"=>[[wire.(S) for S in layer] for layer in evens],
        "boundaries"=>boundaries,"generators"=>wire.(C.generators),
        "name"=>"$(C.group)($(P.rank),O), D=$(P.discriminant), native Julia",
        "discriminant"=>P.discriminant,"rank"=>P.rank,"groupKind"=>string(C.group))
    GAP.Globals.HAPArithmeticVoronoiComplex(GAP.GapObj(data;recursive=true),resolution_length)
end

"""A fresh native Julia Voronoi complex as a HAP non-free resolution object."""
function hap_voronoi_complex(D::Int;group=:GL,rank=2,resolution_length=100,progress=false)
    hap_voronoi_complex(voronoi_request(D;group,rank);resolution_length,progress)
end

"""Construct a HAP complex from a field/order/lattice/group request record."""
function hap_voronoi_complex(request::VoronoiRequestRecord;resolution_length=100,progress=false)
    r=normalize_voronoi_request(request)
    group=Symbol(r["group"])
    group==:SL && r["field"]["discriminant"]>0 && throw(ArgumentError("real-quadratic SL complexes are not implemented"))
    resolution_length isa Integer && resolution_length>0 || throw(ArgumentError("resolution_length must be positive"))
    E=enumerate_perfect_forms(arithmetic_problem(r);backend=:native,progress)
    hap_complex(voronoi_complex(E;group,progress);resolution_length)
end

"""
Register the in-process Julia provider with HAP's optional backend dispatcher.
Afterwards use VoronoiComplex(request : backend := "julia") in GAP. The
quadratic GL/SL convenience constructors still accept squarefree radicands.
"""
function enable_hap_julia_backend(;progress=false,resolution_length=100)
    resolution_length isa Integer && resolution_length>0 || throw(ArgumentError("resolution_length must be positive"))
    GAP.Packages.load("hap") || error("HAP is not installed")
    GAP.Globals.Read(GAP.GapObj(joinpath(PACKAGE_ROOT,"gap","VoronoiBackends.g")))
    provider=function(request)
        hap_voronoi_complex(GAP.gap_to_julia(request;recursive=true);progress,resolution_length)
    end
    GAP.Globals.RegisterVoronoiBackend(GAP.GapObj("julia"),GAP.GapObj(provider))
    nothing
end
