module ArithmeticVoronoi

using Oscar, JSON3, SHA, Dates
import Hecke

export ArithmeticVoronoiProblem, arithmetic_problem, classical_problem,
       trace_gram, form_from_trace, encode_matrix, decode_matrix, check_group_matrix,
       field_vector, coordinate_vector, evaluation_rays, positive_definite,
       shortest_data, equivalence, finite_stabilizer, backend_gram, internal_gram,
       internal_matrix, internal_vector, PerfectFormRecord, NeighbourRecord,
       ArithmeticVoronoiEnumeration, GroupConstraintError

const SCHEMA="arithmetic_voronoi_v1"
const ACTION="columns; F^g=transpose(g)*F*g; (F^g)^h=F^(g*h)"
const HERMITIAN_ACTION="columns; F^g=conjugate_transpose(g)*F*g; (F^g)^h=F^(g*h)"

struct ArithmeticVoronoiProblem
    field::Any
    order::Any
    generator::Any
    degree::Int
    rank::Int
    discriminant::Int
    root_sum::Int
    root_product::Int
    multiplication::QQMatrix
    form_basis::Vector{QQMatrix}
    action_convention::String
    group_constraints::String
end

struct GroupConstraintError <: Exception
    reason::String
end
Base.showerror(io::IO,e::GroupConstraintError)=print(io,e.reason)

is_hermitian(P::ArithmeticVoronoiProblem)=P.discriminant<0
form_conjugate(P,x)=is_hermitian(P) ? P.field(coeff(x,0))+coeff(x,1)*(P.root_sum-P.generator) : x
form_adjoint(P,A)=is_hermitian(P) ? transpose(map_entries(x->form_conjugate(P,x),A)) : transpose(A)
act_form(P,F,g)=form_adjoint(P,g)*F*g
# Historical symmetric-form helper; arithmetic code uses the explicit problem.
act_form(F,g)=transpose(g)*F*g
vector_outer(P,v)=matrix(P.field,[v[i]*form_conjugate(P,v[j]) for i in 1:P.rank,j in 1:P.rank])
form_lattice(P,F)=is_hermitian(P) ? hermitian_lattice(P.field;gram=F) : quadratic_lattice(P.field;gram=F)

function trace_block(s,p,a0,a1)
    t0,t1,t2,t3=2,s,s*s-2p,s*s*s-3s*p
    matrix(QQ,[a0*t0+a1*t1 a0*t1+a1*t2;a0*t1+a1*t2 a0*t2+a1*t3])
end

function trace_gram(P::ArithmeticVoronoiProblem,F)
    size(F)==(P.rank,P.rank) && F==form_adjoint(P,F) || throw(ArgumentError("form must be symmetric/Hermitian of the specified rank"))
    P.degree==1 && return change_base_ring(QQ,F)
    A=zero_matrix(QQ,2P.rank,2P.rank)
    if is_hermitian(P)
        basis=(P.field(1),P.generator)
        for i in 1:P.rank,j in 1:P.rank,a in 1:2,b in 1:2
            A[2i-2+a,2j-2+b]=tr(form_conjugate(P,basis[a])*F[i,j]*basis[b])
        end
        return A
    end
    for i in 1:P.rank,j in 1:P.rank
        a=P.field(F[i,j]);block=trace_block(P.root_sum,P.root_product,coeff(a,0),coeff(a,1))
        A[2i-1:2i,2j-1:2j]=block
    end
    A
end

function make_problem(K,O,w,d,D,s,p,r=2)
    r isa Integer && !(r isa Bool) && 2<=r<=typemax(Int) || throw(ArgumentError("rank must be an integer at least 2"))
    r=Int(r)
    J=d==1 ? zero_matrix(QQ,1,1) : matrix(QQ,[0 -p;1 s])
    E=block_diagonal_matrix([J for _ in 1:r])
    action=D<0 ? HERMITIAN_ACTION : ACTION
    placeholder=ArithmeticVoronoiProblem(K,O,w,d,r,D,s,p,E,QQMatrix[],action,"GL($r,O) acting K-linearly on the free lattice O^$r; retain -I")
    basis=QQMatrix[]
    for i in 1:r,j in i:r,k in 0:d-1
        D<0 && i==j && k==1 && continue
        F=zero_matrix(K,r,r);a=k==0 ? K(1) : w
        F[i,j]=a;F[j,i]=form_conjugate(placeholder,a)
        push!(basis,trace_gram(placeholder,F))
    end
    ArithmeticVoronoiProblem(K,O,w,d,r,D,s,p,E,basis,action,placeholder.group_constraints)
end

"""K-linear GL(rank,O) forms for a fundamental discriminant D; Hermitian when D<0."""
function arithmetic_problem(D::Int;rank=2)
    (D<0 || D>1) && Hecke.is_fundamental_discriminant(D) || throw(ArgumentError("a fundamental quadratic discriminant is required"))
    s=mod(D,4)
    p=div(s*s-D,4)
    R,x=polynomial_ring(QQ,"t")
    K,w=number_field(x^2-s*x+p,"w";cached=false)
    O=maximal_order(K)
    discriminant(O)==D && w in O || error("generator does not give the requested maximal order")
    # The discriminant of (1,w) equals that of O, so Z[w]=O exactly.
    s*s-4p==D || error("integral basis discriminant mismatch")
    make_problem(K,O,w,2,D,s,p,rank)
end
classical_problem(;rank=2)=make_problem(QQ,ZZ,QQ(0),1,1,0,0,rank)

function form_from_trace(P::ArithmeticVoronoiProblem,A)
    A=change_base_ring(QQ,A)
    size(A)==(P.rank*P.degree,P.rank*P.degree) && A==transpose(A) || throw(ArgumentError("invalid trace Gram dimensions or symmetry"))
    P.degree==1 && return A
    F=zero_matrix(P.field,P.rank,P.rank)
    t1=P.root_sum;t2=t1*t1-2P.root_product;delta=2t2-t1*t1
    for i in 1:P.rank,j in 1:P.rank
        b0,b1=A[2i-1,2j-1],A[2i-1,2j]
        a0=(t2*b0-t1*b1)/delta;a1=(-t1*b0+2b1)/delta
        F[i,j]=P.field(a0)+a1*P.generator
    end
    trace_gram(P,F)==A || throw(ArgumentError("matrix lies outside the specified trace-form space"))
    F
end

function encode_matrix(P::ArithmeticVoronoiProblem,g)
    size(g)==(P.rank,P.rank) || throw(ArgumentError("matrix rank mismatch"))
    P.degree==1 && return change_base_ring(QQ,g)
    J=P.multiplication[1:2,1:2];I2=identity_matrix(QQ,2)
    A=zero_matrix(QQ,2P.rank,2P.rank)
    for i in 1:P.rank,j in 1:P.rank
        a=P.field(g[i,j])
        A[2i-1:2i,2j-1:2j]=coeff(a,0)*I2+coeff(a,1)*J
    end
    A
end

function decode_matrix(P::ArithmeticVoronoiProblem,A;require_group=true)
    A=change_base_ring(QQ,A)
    size(A)==(P.rank*P.degree,P.rank*P.degree) || throw(GroupConstraintError("matrix rank mismatch"))
    if P.degree==1
        g=A
    else
        A*P.multiplication==P.multiplication*A || throw(GroupConstraintError("matrix is not K-linear: it does not commute with multiplication by w"))
        g=zero_matrix(P.field,P.rank,P.rank)
        for i in 1:P.rank,j in 1:P.rank
            g[i,j]=P.field(A[2i-1,2j-1])+A[2i,2j-1]*P.generator
        end
        encode_matrix(P,g)==A || throw(GroupConstraintError("field block reconstruction failed"))
    end
    if require_group
        all(isone(denominator(a)) for a in A) || throw(GroupConstraintError("matrix does not preserve the integral lattice"))
        abs(det(A))==1 || throw(GroupConstraintError("matrix is not unimodular over Z"))
        determinant=det(g)
        unit=P.degree==1 ? abs(determinant)==1 : determinant in P.order && abs(norm(determinant))==1
        unit || throw(GroupConstraintError("determinant is not a unit of the specified order"))
    end
    g
end
check_group_matrix(P,A)=(decode_matrix(P,A);true)

function field_vector(P,z)
    length(z)==P.rank*P.degree || throw(ArgumentError("vector dimension mismatch"))
    P.degree==1 && return [P.field(a) for a in z]
    [P.field(z[2i-1])+z[2i]*P.generator for i in 1:P.rank]
end
function coordinate_vector(P,v)
    length(v)==P.rank || throw(ArgumentError("vector rank mismatch"))
    P.degree==1 && return [QQ(a) for a in v]
    [QQ(coeff(P.field(v[i]),k)) for i in 1:P.rank for k in 0:1]
end

function positive_definite(A)
    A==transpose(A) && all(det(A[1:i,1:i])>0 for i in 1:nrows(A))
end
quadratic_value(A,v)=sum(QQ(v[i])*A[i,j]*QQ(v[j]) for i in 1:nrows(A),j in 1:ncols(A))
function evaluation_rays(P,vectors)
    matrix(QQ,[quadratic_value(B,v) for v in vectors,B in P.form_basis])
end
function signed_key(v)
    a=[BigInt(ZZ(x)) for x in v];i=findfirst(!iszero,a)
    isnothing(i) && throw(ArgumentError("zero is not a minimal vector"))
    a[i]<0 && (a=-a)
    Tuple(a)
end
function shortest_data(A)
    positive_definite(A) || throw(ArgumentError("positive-definite Gram matrix required"))
    L=integer_lattice(;gram=change_base_ring(QQ,A))
    V=shortest_vectors(L);m=minimum(L)
    keys=sort!(unique(signed_key(v) for v in V))
    vectors=[collect(v) for v in keys]
    all(quadratic_value(A,v)==m for v in vectors) || error("short-vector coordinate/length convention mismatch")
    (;minimum=m,vectors,lattice=L,sign_convention="one representative of +/-v, first nonzero coordinate positive; squared trace lengths")
end

function equivalence(P,F,H)
    L=form_lattice(P,H);M=form_lattice(P,F)
    ok,h=Hecke.is_isometric_with_isometry(L,M)
    ok || return nothing
    g=form_adjoint(P,h)
    check_group_matrix(P,encode_matrix(P,g))
    act_form(P,F,g)==H || error("Hecke isometry action convention mismatch")
    g
end
function finite_stabilizer(P,F)
    L=form_lattice(P,F)
    gens=[form_adjoint(P,g) for g in Hecke.automorphism_group_generators(L)]
    central=-identity_matrix(P.field,P.rank)
    central in gens || push!(gens,central)
    for g in gens
        check_group_matrix(P,encode_matrix(P,g))
        act_form(P,F,g)==F || error("Hecke stabilizer action convention mismatch")
    end
    (;generators=gens,order=BigInt(Hecke.automorphism_group_order(L)),lattice=L)
end

backend_permutation(P)=[P.degree*(i-1)+k for k in 1:P.degree for i in 1:P.rank]
backend_gram(P,A)=A[backend_permutation(P),backend_permutation(P)]
internal_gram(P,A)=A[invperm(backend_permutation(P)),invperm(backend_permutation(P))]
internal_matrix(P,A)=decode_matrix(P,internal_gram(P,transpose(A)))
internal_vector(P,v)=[BigInt(ZZ(x)) for x in v[invperm(backend_permutation(P))]]

struct PerfectFormRecord
    form::Any
    trace_gram::QQMatrix
    minimum::QQFieldElem
    minimal_vectors::Vector{Vector{BigInt}}
    evaluation_rays::QQMatrix
    stabilizer::Vector{Any}
    stabilizer_order::BigInt
    facets::Vector{Any}
end
struct NeighbourRecord
    source::Int
    facet::Vector{Int}
    target::Int
    neighbour_gram::QQMatrix
    transporter::Any
    facet_map::Vector{Int}
end
mutable struct ArithmeticVoronoiEnumeration
    problem::ArithmeticVoronoiProblem
    representatives::Vector{PerfectFormRecord}
    neighbours::Vector{NeighbourRecord}
    pending_facets::Vector{Any}
    completion_status::Symbol
    provenance::Dict{String,Any}
end

include("Validation.jl")
include("Backend.jl")
include("Requests.jl")
include("TwoCells.jl")
include("NativeBackend.jl")
include("Presentations.jl")
include("SLExtraction.jl")
include("Simplification.jl")
include("Fields.jl")
include("Projective.jl")
include("FieldPipeline.jl")
include("GLPipeline.jl")
include("Words.jl")
include("Dataset.jl")
include("HAP.jl")

end
