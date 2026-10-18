export BackendConfig, BackendError, default_backend, run_exporter, cache_key,
       enumerate_perfect_forms, load_checkpoint, write_enumeration, enumeration_data

const BACKEND_COMMIT="cfcd79dbddf67a494a1c71b776174b190da5209a"
const PACKAGE_ROOT=normpath(joinpath(@__DIR__,".."))
struct BackendError <: Exception
    kind::Symbol
    message::String
    run_directory::String
end
BackendError(kind,message)=BackendError(kind,message,"")
Base.showerror(io::IO,e::BackendError)=print(io,e.kind,": ",e.message,isempty(e.run_directory) ? "" : " (diagnostics: $(e.run_directory))")

file_sha(path)=bytes2hex(open(SHA.sha256,path))
struct BackendConfig
    executable::String
    executable_sha256::String
    commit::String
end
function BackendConfig(path::AbstractString;expected_sha256=nothing,commit=BACKEND_COMMIT)
    isabspath(path) || throw(BackendError(:invalid_executable,"configure one absolute executable path"))
    isfile(path) && !isnothing(Sys.which(path)) || throw(BackendError(:missing_executable,"missing or nonexecutable backend: $path"))
    path=realpath(path);actual=file_sha(path)
    isnothing(expected_sha256) || actual==expected_sha256 || throw(BackendError(:incompatible_backend,"executable hash differs from its configured build"))
    BackendConfig(path,actual,commit)
end
function default_backend()
    path=joinpath(PACKAGE_ROOT,"backend/config.json")
    isfile(path) || throw(BackendError(:missing_executable,"build the exporter and create backend/config.json, or pass BackendConfig(absolute_path)"))
    data=JSON3.read(read(path,String),Dict{String,Any})
    BackendConfig(data["executable"];expected_sha256=data["executable_sha256"],commit=data["commit"])
end
function default_backend(P::ArithmeticVoronoiProblem)
    P.discriminant in (12,13,40) || return default_backend()
    path=joinpath(PACKAGE_ROOT,"backend/config_fields.json")
    isfile(path) || throw(BackendError(:missing_executable,"build and configure the extended-field exporter"))
    data=JSON3.read(read(path,String),Dict{String,Any})
    BackendConfig(data["executable"];expected_sha256=data["executable_sha256"],commit=data["commit"])
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

function namelist(P;arithmetic="gmp")
    arithmetic in ("gmp","flint","unsupported") || throw(ArgumentError("arithmetic name contains unsupported input characters"))
    """
    &TSPACE
      TypeTspace = "$(P.degree==1 ? "Classic" : "RealQuad")"
      ClassicDim = 2
      RealImagDim = 2
      RealImagSum = $(P.root_sum)
      RealImagProd = $(P.root_product)
      SuperMatMethod = "NotNeeded"
      ListComm = "$(P.degree==1 ? "Trivial" : "Use_realimag")"
      PtGroupMethod = "Trivial"
      FileListSubspaces = "unset"
    /
    &DATA
      arithmetic = "$arithmetic"
      OnlyWellRounded = T
      ComputeComplex = F
      ComputeBoundary = F
      ComputePreImage = F
      CacheFile = "none"
    /
    """
end

function child_error(stderr)
    allowed=Set(["unsupported_arithmetic","unsupported_operation","malformed_input","unwritable_output","backend_failure"])
    for line in reverse(split(stderr,'\n'))
        startswith(line,"AV_ERROR ") || continue
        try
            data=JSON3.read(line[10:end],Dict{String,Any})
            if data["schema"]=="av_error_v1" && data["code"] in allowed
                return Symbol(data["code"]),string(data["message"])
            end
        catch
        end
    end
    :child_exit,"backend exited without a valid structured error record"
end

function resident_bytes(pid,ps)
    text=read(pipeline(ignorestatus(`$ps -o rss= -p $pid`);stderr=devnull),String)
    rss=tryparse(Int,strip(text))
    isnothing(rss) ? nothing : 1024rss
end

function stop_backend(process)
    process_running(process) || return
    kill(process,Base.SIGTERM)
    for _ in 1:20
        process_exited(process) && return
        sleep(.01)
    end
    process_running(process) && kill(process,Base.SIGKILL)
end

"""Run an executable with private files, an outer deadline and bounded memory.

The default memory policy is a hard address-space limit on Linux and a sampled
resident-memory watchdog on macOS. The resident policy can also be selected on
Linux, including for ASan executables whose virtual shadow memory is large.
"""
function run_exporter(P,B::BackendConfig;work_dir=joinpath(PACKAGE_ROOT,"results/runs"),timeout=60.0,
                      input_text=namelist(P),output_path=nothing,memory_limit=10*1024^3,
                      memory_policy=:auto)
    isfinite(timeout) && timeout>0 || throw(ArgumentError("positive finite timeout required"))
    memory_limit isa Integer && 0<memory_limit<=typemax(Int) || throw(ArgumentError("positive memory limit in bytes required"))
    memory_policy in (:auto,:address_space,:resident) || throw(ArgumentError("memory policy must be :auto, :address_space or :resident"))
    work_dir=abspath(work_dir)
    try mkpath(work_dir) catch err
        throw(BackendError(:unwritable_output,sprint(showerror,err)))
    end
    folder=mktempdir(work_dir;prefix="request ",cleanup=false)
    input=joinpath(folder,"input.nml");write(input,input_text)
    token=bytes2hex(sha256(input_text*string(time_ns(),getpid())))
    output=isnothing(output_path) ? joinpath(folder,"output.json") : abspath(output_path)
    stdout_path=joinpath(folder,"stdout.log");stderr_path=joinpath(folder,"stderr.log")
    report=Dict{String,Any}("status"=>"prepared","completed"=>false,"request_token"=>token,
        "input_sha256"=>file_sha(input),"backend_commit"=>B.commit,"executable"=>B.executable,
        "executable_sha256"=>B.executable_sha256,"timeout_seconds"=>timeout,"memory_limit_bytes"=>memory_limit,
        "output"=>output,"started_utc"=>string(now(UTC)))
    report_path=joinpath(folder,"process.json")
    function fail(kind,message)
        report["status"]=string(kind);report["error"]=message
        atomic_json(report_path,report)
        throw(BackendError(kind,message,folder))
    end
    isfile(B.executable) || fail(:missing_executable,"configured executable disappeared")
    file_sha(B.executable)==B.executable_sha256 || fail(:incompatible_backend,"configured executable changed")
    ispath(output) && fail(:stale_output,"output path already exists before child launch")
    try
        probe,probe_io=mktemp(dirname(output);cleanup=false);close(probe_io);rm(probe)
    catch err
        fail(:unwritable_output,sprint(showerror,err))
    end
    Sys.islinux() || Sys.isapple() || fail(:unsupported_platform,"backend process control supports Linux and macOS")
    policy=memory_policy==:auto ? (Sys.islinux() ? :address_space : :resident) : memory_policy
    argv=[B.executable,input,output,token]
    ps=nothing
    if policy==:address_space
        Sys.islinux() || fail(:unsupported_platform,"address-space limiting requires Linux; use :resident on macOS")
        limiter=Sys.which("prlimit")
        isnothing(limiter) && fail(:unsupported_platform,"install util-linux for prlimit, or use memory_policy=:resident")
        argv=[limiter,"--as=$memory_limit","--",argv...]
    else
        ps=Sys.which("ps")
        isnothing(ps) && fail(:unsupported_platform,"resident-memory monitoring requires ps")
    end
    report["memory_policy"]=string(policy)
    report["maximum_sampled_rss_bytes"]=0
    report["argv"]=argv;atomic_json(report_path,report)
    started=time_ns();failure=nothing;code=0;signal=0
    try
        open(stdout_path,"w") do stdout_io
            open(stderr_path,"w") do stderr_io
                process=run(pipeline(Cmd(Cmd(argv);dir=folder);stdout=stdout_io,stderr=stderr_io);wait=false)
                try
                    while process_running(process)
                        if (time_ns()-started)/1e9>timeout
                            failure=(:timeout,"outer deadline terminated the backend; no completion checkpoint published")
                            break
                        end
                        if !isnothing(ps)
                            rss=resident_bytes(getpid(process),ps)
                            if isnothing(rss)
                                process_running(process) || break
                                failure=(:memory_monitor_error,"could not read backend resident memory; backend stopped")
                                break
                            end
                            report["maximum_sampled_rss_bytes"]=max(report["maximum_sampled_rss_bytes"],rss)
                            if rss>memory_limit
                                failure=(:memory_limit,"resident-memory watchdog terminated the backend; no completion checkpoint published")
                                break
                            end
                        end
                        sleep(min(.02,timeout/10))
                    end
                finally
                    stop_backend(process)
                    wait(process)
                end
                code=process.exitcode;signal=process.termsignal
            end
        end
    catch err
        err isa BackendError && rethrow()
        fail(:launch_error,sprint(showerror,err))
    end
    report["elapsed_seconds"]=(time_ns()-started)/1e9
    report["exitcode"]=code;report["termination_signal"]=signal
    isnothing(failure) || fail(failure...)
    signal!=0 && fail(:terminated,"backend terminated by signal $signal")
    if code!=0
        kind,message=child_error(read(stderr_path,String));fail(kind,message)
    end
    isfile(output) || fail(:missing_output,"successful child exit without the requested output")
    raw=read(output,String)
    data=try JSON3.read(raw,Dict{String,Any}) catch err
        fail(:truncated_output,sprint(showerror,err))
    end
    get(data,"request_token",nothing)==token || fail(:stale_output,"fresh output has the wrong request token")
    get(data,"schema",nothing)=="av_polyhedral_export_v1" || fail(:incompatible_schema,"unsupported exporter schema")
    get(data,"completed",false)===true || fail(:incomplete_output,"backend did not declare completed enumeration")
    get(data,"backend_commit",nothing)==B.commit || fail(:incompatible_backend,"exported backend version does not match configuration")
    report["status"]="exported_pending_math_validation";report["raw_output_sha256"]=bytes2hex(sha256(raw))
    atomic_json(report_path,report)
    (;data,raw,run_directory=folder,process_info=report)
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

function enumerate_perfect_forms(P;backend=:polyhedral_common,config=nothing,
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
