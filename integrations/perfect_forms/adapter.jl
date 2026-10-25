module PerfectFormsAdapter

using ArithmeticVoronoi, JSON3, SHA, Dates
import ArithmeticVoronoi: SCHEMA, problem_descriptor, validate_export, atomic_json
import PerfectForms
import PerfectForms: BackendConfig, BackendError, default_backend, run_exporter

export BackendConfig, BackendError, default_backend, run_exporter, cache_key,
       load_checkpoint

# Load explicitly from the integration environment; the registered engine has
# no dependency on the unregistered C++ package.
default_backend(::ArithmeticVoronoiProblem)=default_backend()
function run_exporter(P::ArithmeticVoronoiProblem,B::BackendConfig;
                      input_text=namelist(P),work_dir=joinpath(pwd(),"results/runs"),kwargs...)
    run_exporter(input_text,B;work_dir,kwargs...)
end

backend_options(P)=(;arithmetic="gmp",ListComm=P.degree==1 ? "Trivial" : "Use_realimag",PtGroupMethod="Trivial (-I)",
    seed=1,computed_levels=["perfect_forms","facet_orbits","labelled_neighbours"],schema="av_polyhedral_export_v1")
function cache_key(P,B::BackendConfig)
    descriptor=(;schema=SCHEMA,problem=problem_descriptor(P),backend_commit=B.commit,
                executable_sha256=B.executable_sha256,options=backend_options(P))
    bytes2hex(sha256(JSON3.write(descriptor)))
end
function namelist(P;arithmetic="gmp")
    P.rank==2 || throw(ArgumentError("the C++ request adapter supports rank 2; use backend=:native for higher ranks"))
    PerfectForms.perfect_form_request(P.discriminant;arithmetic)
end

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

function ArithmeticVoronoi.enumerate_external(P::ArithmeticVoronoiProblem,::Val{:polyhedral_common};config=nothing,
        cache_dir=joinpath(pwd(),"results/cache"),work_dir=joinpath(pwd(),"results/runs"),
        timeout=60.0,force=true,progress=false)
    P.rank==2 || throw(ArgumentError("the C++ request adapter supports rank 2; use backend=:native for higher ranks"))
    B=isnothing(config) ? default_backend(P) : config
    key=cache_key(P,B);path=joinpath(cache_dir,key*".json")
    if isfile(path) && !force
        @info "Loading a validated perfect-form checkpoint" path
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

end # module
