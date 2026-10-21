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
    isempty(V) && return nothing
    u=field_vector(P,V[1])
    for j in 2:length(V)
        v=field_vector(P,V[j])
        u[1]*v[2]!=u[2]*v[1] && return matrix(P.field,[u[1] v[1];u[2] v[2]])
    end
    nothing
end

function make_voronoi_cell(P,V,F)
    vectors=[collect(v) for v in cell_vector_key(V)]
    ensure(!isnothing(cell_basis(P,vectors)),"cell minimum set does not span K^2")
    rays=evaluation_rays(P,vectors)
    _,T=nullspace(rays)
    S=sum(vector_outer(P,field_vector(P,v)) for v in vectors)
    ensure(positive_definite(trace_gram(P,S)),"cell has no positive interior in the dual cone")
    VoronoiCell(vectors,F,T,inv(S))
end

function transform_cell(P,c,g)
    make_voronoi_cell(P,cell_image(P,g,c.vectors),act_form(P,c.point,inv(g)))
end

# Exhaustive rank-two transporter search: an integral linear map is determined
# by the images of two K-independent vectors. Both signs of every target vector
# are used; this retains all scalar units, including central -I.
function cell_maps(P,source,target; determinant_one=false, first_only=false)
    B=cell_basis(P,source)
    ensure(!isnothing(B),"transporter source is not well rounded")
    inverse=inv(B)
    W=Set(cell_vector_key(target))
    signed=[field_vector(P,v) for v in [target;[-v for v in target]]]
    result=Any[]
    for u in signed,v in signed
        iszero(u[1]*v[2]-u[2]*v[1]) && continue
        g=matrix(P.field,[u[1] v[1];u[2] v[2]])*inverse
        delta=det(g)
        determinant_one && delta!=1 && continue
        (P.degree==1 ? abs(delta)==1 : abs(norm(delta))==1) || continue
        M=encode_matrix(P,g)
        all(x->isone(denominator(x)),M) || continue
        all(z->signed_key(apply_integral(M,z)) in W,source) || continue
        push!(result,g)
        first_only && return result
    end
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
equivariant boundaries from a complete rank-two perfect-form enumeration.
Supports quadratic GL and imaginary-quadratic SL on the free lattice O^2.
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
        reps=[diagonal_matrix(P.field,[units.torsion_generator^j,one(P.field)]) for j in 0:units.torsion_order-1]
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
    for (k,layer) in enumerate(C.cells),(i,c) in enumerate(layer)
        group=C.stabilizers[k][i]
        ensure(-identity_matrix(P.field,2) in group,"cell stabilizer lost central -I")
        for g in group
            check_group_matrix(P,encode_matrix(P,g))
            ensure(cell_vector_key(cell_image(P,g,c.vectors))==cell_vector_key(c.vectors),"invalid cell stabilizer")
            cell_orientation(P,c,g)
        end
        for (j,g) in C.boundaries[k][i]
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
    # All these groups contain -I. Reuse the verified finite-group reduction.
    finite_group_record(P,elements).generators
end

function hap_complex(C::VoronoiComplex;resolution_length=100)
    resolution_length isa Integer && resolution_length>0 || throw(ArgumentError("resolution_length must be positive"))
    validate_voronoi_complex(C)
    GAP.Packages.load("hap") || error("HAP is not installed; use GAP.Packages.install(\"hap\")")
    GAP.Globals.Read(GAP.GapObj(joinpath(PACKAGE_ROOT,"gap","ArithmeticVoronoi.g")))
    P=C.enumeration.problem
    elts=Any[identity_matrix(P.field,2)]
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
        "name"=>"$(C.group)(2,O), D=$(P.discriminant), native Julia",
        "discriminant"=>P.discriminant,"groupKind"=>string(C.group))
    GAP.Globals.HAPArithmeticVoronoiComplex(GAP.GapObj(data;recursive=true),resolution_length)
end

"""A fresh native Julia Voronoi complex as a HAP non-free resolution object."""
function hap_voronoi_complex(D::Int;group=:GL,resolution_length=100,progress=false)
    group in (:GL,:SL) || throw(ArgumentError("choose group=:GL or :SL"))
    group==:SL && D>0 && throw(ArgumentError("real-quadratic SL complexes are not implemented"))
    E=enumerate_perfect_forms(arithmetic_problem(D);backend=:native,progress)
    hap_complex(voronoi_complex(E;group,progress);resolution_length)
end

"""
Register the in-process Julia provider with HAP's optional backend dispatcher.
Afterwards use VoronoiComplexGL(2,d : backend := "julia") in GAP. HAP takes
squarefree radicands d; hap_voronoi_complex takes fundamental discriminants D.
"""
function enable_hap_julia_backend(;progress=false)
    GAP.Packages.load("hap") || error("HAP is not installed")
    GAP.Globals.Read(GAP.GapObj(joinpath(PACKAGE_ROOT,"gap","VoronoiBackends.g")))
    provider=function(kind,n,d)
        n=Int(n);d=Int(d);kind=Symbol(String(GAP.gap_to_julia(kind)))
        n==2 || throw(ArgumentError("the Julia Voronoi backend supports rank 2"))
        d!=0 && d!=1 && is_squarefree(ZZ(d)) || throw(ArgumentError("HAP requires a squarefree quadratic radicand"))
        D=mod(d,4)==1 ? d : 4d
        hap_voronoi_complex(D;group=kind,progress)
    end
    GAP.Globals.RegisterVoronoiBackend(GAP.GapObj("julia"),GAP.GapObj(provider))
    nothing
end
