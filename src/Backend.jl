export BackendConfig, BackendError, default_backend, run_exporter, cache_key,
       enumerate_perfect_forms, load_checkpoint, write_enumeration, enumeration_data

import PerfectForms
import PerfectForms: BackendConfig, BackendError, default_backend, run_exporter,
                     enumerate_perfect_forms, BACKEND_COMMIT, file_sha, atomic_json

const PACKAGE_ROOT=normpath(joinpath(@__DIR__,".."))

# Preserve the problem-based API; executable discovery and process control belong
# to PerfectForms, which has no dependency on our arithmetic types.
default_backend(::ArithmeticVoronoiProblem)=default_backend()
function run_exporter(P::ArithmeticVoronoiProblem,B::BackendConfig;
                      input_text=namelist(P),work_dir=joinpath(PACKAGE_ROOT,"results/runs"),kwargs...)
    run_exporter(input_text,B;work_dir,kwargs...)
end

fraction_string(q)=begin
    q=QQ(q);d=denominator(q);isone(d) ? string(numerator(q)) : string(numerator(q),"/",d)
end
matrix_strings(A)=[[fraction_string(A[i,j]) for j in 1:ncols(A)] for i in 1:nrows(A)]
field_matrix_data(P,A)=[[[fraction_string(P.degree==1 ? A[i,j] : coeff(A[i,j],k)) for k in 0:P.degree-1] for j in 1:ncols(A)] for i in 1:nrows(A)]
function problem_descriptor(P)
    (;discriminant=P.discriminant,degree=P.degree,rank=P.rank,
      polynomial_coefficients=P.degree==1 ? [0,1] : [P.root_product,-P.root_sum,1],
      order="maximal",order_basis=P.degree==1 ? ["1"] : ["1","w"],
      lattice="free O^2",coordinate_order=P.degree==1 ? ["e1","e2"] : ["e1","w*e1","e2","w*e2"],
      multiplication=matrix_strings(P.multiplication),form_basis=matrix_strings.(P.form_basis),
      group=P.group_constraints,action=P.action_convention)
end
backend_options(P)=(;arithmetic="gmp",ListComm=P.degree==1 ? "Trivial" : "Use_realimag",PtGroupMethod="Trivial (-I)",
    seed=1,computed_levels=["perfect_forms","facet_orbits","labelled_neighbours"],schema="av_polyhedral_export_v1")
function cache_key(P,B::BackendConfig)
    descriptor=(;schema=SCHEMA,problem=problem_descriptor(P),backend_commit=B.commit,
                executable_sha256=B.executable_sha256,options=backend_options(P))
    bytes2hex(sha256(JSON3.write(descriptor)))
end
namelist(P;arithmetic="gmp")=PerfectForms.perfect_form_request(P.discriminant;arithmetic)

function load_checkpoint(P,B,path)
    checkpoint=try JSON3.read(read(path,String),Dict{String,Any}) catch err
        throw(BackendError(:truncated_cache,sprint(showerror,err)))
    end
    get(checkpoint,"schema",nothing)=="av_checkpoint_v1" || throw(BackendError(:incompatible_schema,"unsupported cache version"))
    get(checkpoint,"completed",false)===true || throw(BackendError(:incomplete_cache,"cache was not atomically committed after validation"))
    required=("problem","backend_commit","executable_sha256","options","cache_key","raw_export_json","raw_export_sha256")
    all(haskey(checkpoint,k) for k in required) && checkpoint["problem"] isa AbstractDict && checkpoint["raw_export_json"] isa AbstractString ||
        throw(BackendError(:incompatible_cache,"missing or malformed cache envelope fields"))
    expected=JSON3.read(JSON3.write(problem_descriptor(P)),Dict{String,Any})
    problem=checkpoint["problem"]
    get(problem,"discriminant",nothing)==P.discriminant && get(problem,"degree",nothing)==P.degree || throw(BackendError(:wrong_field,"cache belongs to a different number field"))
    get(problem,"group",nothing)==P.group_constraints && get(problem,"action",nothing)==P.action_convention || throw(BackendError(:wrong_subgroup,"cache group/action constraints differ"))
    problem==expected || throw(BackendError(:incompatible_problem,"order, basis or lattice data differ"))
    checkpoint["backend_commit"]==B.commit && checkpoint["executable_sha256"]==B.executable_sha256 || throw(BackendError(:incompatible_backend,"cache comes from a different backend build"))
    checkpoint["options"]==JSON3.read(JSON3.write(backend_options(P)),Dict{String,Any}) || throw(BackendError(:incompatible_options,"cache options or commuting constraints differ"))
    checkpoint["cache_key"]==cache_key(P,B) || throw(BackendError(:incompatible_cache,"cache key mismatch"))
    raw=checkpoint["raw_export_json"]
    bytes2hex(sha256(raw))==checkpoint["raw_export_sha256"] || throw(BackendError(:corrupt_cache,"export checksum mismatch"))
    data=try JSON3.read(raw,Dict{String,Any}) catch err
        throw(BackendError(:truncated_cache,sprint(showerror,err)))
    end
    get(data,"backend_commit",nothing)==B.commit || throw(BackendError(:incompatible_backend,"raw cached output declares a different backend commit"))
    start=time_ns()
    E=try validate_export(P,data) catch err
        throw(BackendError(:arithmetic_validation,sprint(showerror,err)))
    end
    E.provenance["from_cache"]=true;E.provenance["cache_key"]=cache_key(P,B)
    E.provenance["executable_sha256"]=B.executable_sha256
    E.provenance["validation_seconds"]=(time_ns()-start)/1e9
    E
end

function enumerate_perfect_forms(P::ArithmeticVoronoiProblem;backend=:polyhedral_common,config=nothing,
        cache_dir=joinpath(PACKAGE_ROOT,"results/cache"),work_dir=joinpath(PACKAGE_ROOT,"results/runs"),
        timeout=60.0,force=false)
    backend==:polyhedral_common || throw(BackendError(:unsupported_backend,"native backend is not implemented; the mathematical records are backend-independent"))
    B=isnothing(config) ? default_backend(P) : config
    key=cache_key(P,B);path=joinpath(cache_dir,key*".json")
    if isfile(path) && !force
        return load_checkpoint(P,B,path)
    end
    run=run_exporter(P,B;work_dir,timeout)
    start=time_ns()
    E=try validate_export(P,run.data) catch err
        report=run.process_info;report["status"]="arithmetic_validation_failed";report["error"]=sprint(showerror,err)
        atomic_json(joinpath(run.run_directory,"process.json"),report)
        throw(BackendError(:arithmetic_validation,sprint(showerror,err),run.run_directory))
    end
    E.provenance["from_cache"]=false;E.provenance["cache_key"]=key
    E.provenance["executable_sha256"]=B.executable_sha256
    E.provenance["run_directory"]=run.run_directory
    E.provenance["backend_seconds"]=run.process_info["elapsed_seconds"]
    E.provenance["validation_seconds"]=(time_ns()-start)/1e9
    checkpoint=(;schema="av_checkpoint_v1",completed=true,cache_key=key,problem=problem_descriptor(P),
        backend_commit=B.commit,executable_sha256=B.executable_sha256,options=backend_options(P),
        raw_export_sha256=bytes2hex(sha256(run.raw)),raw_export_json=run.raw,
        validated_utc=string(now(UTC)),validation=E.provenance)
    atomic_json(path,checkpoint)
    report=run.process_info;report["completed"]=true;report["status"]="validated_complete"
    report["cache_path"]=path;atomic_json(joinpath(run.run_directory,"process.json"),report)
    E
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
      transporter_convention="H=transpose(g)*F_target*g; g maps neighbour minimal vectors to target minimal vectors",
      word_convention="matrix products in written order; negative stabilizer index denotes inverse",
      forms,neighbours,pending_facets=E.pending_facets,provenance=E.provenance)
end
write_enumeration(path,E)=atomic_json(path,enumeration_data(E))
