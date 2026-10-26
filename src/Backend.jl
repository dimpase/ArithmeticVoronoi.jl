export enumerate_perfect_forms, write_enumeration, enumeration_data

const PACKAGE_ROOT=normpath(joinpath(@__DIR__,".."))

file_sha(path)=bytes2hex(open(SHA.sha256,path))
function atomic_json(path,data)
    mkpath(dirname(path));temporary,io=mktemp(dirname(path);cleanup=false)
    try
        JSON3.pretty(io,data);println(io);flush(io);close(io)
        Base.Filesystem.rename(temporary,path)
    catch
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
        rethrow()
    end
    path
end

fraction_string(q)=begin
    q=QQ(q);d=denominator(q);isone(d) ? string(numerator(q)) : string(numerator(q),"/",d)
end
field_scalar_data(P,a)=[fraction_string(P.degree==1 ? a : coeff(a,k)) for k in 0:P.degree-1]
matrix_strings(A)=[[fraction_string(A[i,j]) for j in 1:ncols(A)] for i in 1:nrows(A)]
field_matrix_data(P,A)=[[[fraction_string(P.degree==1 ? A[i,j] : coeff(A[i,j],k)) for k in 0:P.degree-1] for j in 1:ncols(A)] for i in 1:nrows(A)]
function problem_descriptor(P)
    (;discriminant=P.discriminant,degree=P.degree,rank=P.rank,
      polynomial_coefficients=P.degree==1 ? [0,1] : [P.root_product,-P.root_sum,1],
      order="maximal",order_basis=P.degree==1 ? ["1"] : ["1","w"],
      lattice="free O^$(P.rank)",coordinate_order=P.degree==1 ? ["e$i" for i in 1:P.rank] :
        [label for i in 1:P.rank for label in ("e$i","w*e$i")],
      multiplication=matrix_strings(P.multiplication),form_basis=matrix_strings.(P.form_basis),
      group=P.group_constraints,action=P.action_convention)
end
"""
    enumerate_perfect_forms(P::ArithmeticVoronoiProblem; backend=:native, progress=false, ...)

Compute perfect forms, complete K-linear stabilizers and labelled neighbours
from scratch using Julia/Hecke/OSCAR. `progress=true` prints seed and facet
traversal progress. Native enumeration has no deadline or disk cache.

The optional C++ comparison backend requires explicitly loading
`integrations/perfect_forms/adapter.jl` from its separate Julia environment;
then select `backend=:polyhedral_common`. `config`, `cache_dir`, `work_dir`,
`timeout` and `force` apply only to that backend. Its checkpoint loading
requires `force=false`. Both paths return an `ArithmeticVoronoiEnumeration`.
"""
function enumerate_perfect_forms(P::ArithmeticVoronoiProblem;backend=:native,config=nothing,
        cache_dir=joinpath(pwd(),"results/cache"),work_dir=joinpath(pwd(),"results/runs"),
        timeout=60.0,force=true,progress=false)
    backend==:native && return enumerate_native(P;progress)
    enumerate_external(P,Val(backend);config,cache_dir,work_dir,timeout,force,progress)
end

# Optional adapters add methods for their backend tag without changing the
# registered package's dependencies or replacing existing methods.
function enumerate_external(P,backend;kwargs...)
    backend==Val(:polyhedral_common) && throw(ArgumentError(
        "load integrations/perfect_forms/adapter.jl from its Julia environment before selecting backend=:polyhedral_common; see USER_GUIDE.md"))
    throw(ArgumentError("unsupported enumeration backend: $backend"))
end

function enumeration_data(E)
    P=E.problem
    forms=[(;id=i,form_coefficients=field_matrix_data(P,r.form),trace_gram=matrix_strings(r.trace_gram),
        minimum=fraction_string(r.minimum),minimal_vectors=[[string(a) for a in v] for v in r.minimal_vectors],
        minimal_vector_convention="up to sign, first nonzero coordinate positive; signed actions retained",
        evaluation_rays=matrix_strings(r.evaluation_rays),stabilizer_order=string(r.stabilizer_order),
        stabilizer_generators=[field_matrix_data(P,g) for g in r.stabilizer],
        facets=[(;minimal_vector_indices=f.indices,normal=fraction_string.(f.normal),neighbour_orbit=f.neighbour_orbit,
            to_orbit_representative=field_matrix_data(P,f.to_orbit_representative),stabilizer_word=f.stabilizer_word,
            induced_minimal_vector_permutation=f.induced_minimal_vector_permutation) for f in r.facets]) for (i,r) in enumerate(E.representatives)]
    neighbours=[(;source=n.source,source_facet=n.facet,target=n.target,neighbour_trace_gram=matrix_strings(n.neighbour_gram),
        transporter=field_matrix_data(P,n.transporter),signed_target_facet_map=n.facet_map) for n in E.neighbours]
    (;schema=SCHEMA,problem=problem_descriptor(P),completion_status=string(E.completion_status),
      normalization="minimum squared trace length one; one common rational scale",
      transporter_convention=is_hermitian(P) ? "H=conjugate_transpose(g)*F_target*g; g maps neighbour minimal vectors to target minimal vectors" : "H=transpose(g)*F_target*g; g maps neighbour minimal vectors to target minimal vectors",
      word_convention="matrix products in written order; negative stabilizer index denotes inverse",
      forms,neighbours,pending_facets=E.pending_facets,provenance=E.provenance)
end
write_enumeration(path,E)=atomic_json(path,enumeration_data(E))
