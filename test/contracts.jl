using ArithmeticVoronoi, Oscar, JSON3, Test
const CONTRACT_ROOT=mktempdir(joinpath(@__DIR__,"../results");prefix="contract tests ",cleanup=false)
const CONTRACT_AV=ArithmeticVoronoi
const OBSERVED=Dict{String,Any}()

function expect_backend_error(f,label,kind)
    observed=nothing
    try f() catch err;observed=err;end
    @test observed isa BackendError
    if observed isa BackendError
        @test observed.kind==kind
        OBSERVED[label]=Dict("kind"=>string(observed.kind),"expected"=>string(kind),
            "message"=>observed.message,"diagnostics"=>observed.run_directory)
    else
        OBSERVED[label]=Dict("unexpected"=>string(observed))
    end
end

@testset "external boundary, failure behaviour and exact cache contracts" begin
    P=arithmetic_problem(5);B=default_backend()
    @test BackendConfig===CONTRACT_AV.PerfectForms.BackendConfig
    @test run_exporter===CONTRACT_AV.PerfectForms.run_exporter
    @test enumerate_perfect_forms===CONTRACT_AV.PerfectForms.enumerate_perfect_forms
    E=enumerate_perfect_forms(P;config=B)
    @test E.completion_status==:complete
    @test isempty(E.pending_facets)
    @test all(c["unresolved_facets"]==0 for c in E.provenance["closure"])
    @test sort([r.stabilizer_order for r in E.representatives])==[12,20] # verified fixture, not the enumeration stopping rule
    # Public reduction labels and their words/permutations must have the same direction.
    for r in E.representatives,f in r.facets
        word_matrix=identity_matrix(P.field,P.rank)
        for k in f.stabilizer_word
            word_matrix*=k>0 ? r.stabilizer[k] : inv(r.stabilizer[-k])
        end
        @test word_matrix==f.to_orbit_representative
        index=Dict(Tuple(v)=>j for (j,v) in enumerate(r.minimal_vectors))
        R=encode_matrix(P,word_matrix)
        permutation=[index[CONTRACT_AV.signed_key([sum(R[i,j]*v[j] for j in eachindex(v)) for i in eachindex(v)])] for v in r.minimal_vectors]
        @test permutation==f.induced_minimal_vector_permutation
    end
    # The executable itself, the working directory, and the request directory contain spaces.
    spaced=joinpath(CONTRACT_ROOT,"backend directory with spaces");mkpath(spaced)
    copy_path=joinpath(spaced,"exact exporter");cp(B.executable,copy_path);chmod(copy_path,0o700)
    spaced_backend=BackendConfig(copy_path;expected_sha256=B.executable_sha256)
    repeated=enumerate_perfect_forms(P;config=spaced_backend,force=true,
        cache_dir=joinpath(CONTRACT_ROOT,"repeat cache"),work_dir=CONTRACT_ROOT)
    @test repeated.completion_status==:complete
    orbit_map=Int[]
    for record in repeated.representatives
        matching=findall(r->!isnothing(equivalence(P,r.form,record.form)),E.representatives)
        @test length(matching)==1
        push!(orbit_map,only(matching))
    end
    @test sort(orbit_map)==collect(eachindex(E.representatives))
    @test length(repeated.neighbours)==length(E.neighbours)
    OBSERVED["path with spaces and mathematical repeat"]=Dict("status"=>"passed","orbit_map"=>orbit_map)
    # Actual backend interruption, then restart with the same problem/build key.
    restart_cache=joinpath(CONTRACT_ROOT,"killed run cache")
    expect_backend_error("killed-run deadline",:timeout) do
        enumerate_perfect_forms(P;config=B,timeout=.005,cache_dir=restart_cache,work_dir=CONTRACT_ROOT)
    end
    @test !isdir(restart_cache) || isempty(readdir(restart_cache))
    restarted=enumerate_perfect_forms(P;config=B,timeout=60.,cache_dir=restart_cache,work_dir=CONTRACT_ROOT)
    @test restarted.completion_status==:complete && !restarted.provenance["from_cache"]
    resumed=enumerate_perfect_forms(P;config=B,force=false,cache_dir=restart_cache,work_dir=CONTRACT_ROOT)
    @test resumed.provenance["from_cache"] && resumed.completion_status==:complete
    OBSERVED["killed-run restart"]=Dict("status"=>"passed","policy"=>"preserve partial diagnostics; restart enumeration from scratch; only fully validated final caches are reusable")

    cache_path=joinpath(CONTRACT_AV.PACKAGE_ROOT,"results/cache",cache_key(P,B)*".json")
    good=JSON3.read(read(cache_path,String),Dict{String,Any})
    for (label,kind,mutate) in [
        ("wrong field",:wrong_field,d->(d["problem"]["discriminant"]=8)),
        ("wrong subgroup",:wrong_subgroup,d->(d["problem"]["group"]="SL(2,O)")),
        ("wrong order basis",:incompatible_problem,d->(d["problem"]["order_basis"]=["1","sqrt(5)"])),
        ("wrong commuting options",:incompatible_options,d->(d["options"]["ListComm"]="Trivial")),
        ("incompatible schema",:incompatible_schema,d->(d["schema"]="av_checkpoint_v999")),
        ("incomplete cache",:incomplete_cache,d->(d["completed"]=false)),
        ("different backend build",:incompatible_backend,d->(d["executable_sha256"]="0"^64)),
        ("missing cache envelope field",:incompatible_cache,d->delete!(d,"raw_export_json")),
        ("corrupt payload",:corrupt_cache,d->(d["raw_export_json"]*=" "))]
        bad=deepcopy(good);mutate(bad)
        file=joinpath(CONTRACT_ROOT,replace(label," "=>"_")*".json");CONTRACT_AV.atomic_json(file,bad)
        expect_backend_error(label,kind) do;load_checkpoint(P,B,file);end
    end
    truncated=joinpath(CONTRACT_ROOT,"truncated-cache.json");text=read(cache_path,String)
    write(truncated,text[1:div(length(text),2)])
    expect_backend_error("truncated cache",:truncated_cache) do;load_checkpoint(P,B,truncated);end
    # These mutations remain syntactically valid and checksum-consistent: mathematical revalidation must reject them.
    raw=JSON3.read(good["raw_export_json"],Dict{String,Any})
    mutations=[
        ("false boundary assertion",d->(d["boundary_square"]=false)),
        ("incorrect minimum",d->(d["representatives"][1]["minimum"]="2")),
        ("incomplete stabilizer",d->begin
            r=d["representatives"][1];r["stabilizer_generators_row"]=[CONTRACT_AV.matrix_strings(-identity_matrix(QQ,4))];r["stabilizer_order"]="2"
        end),
        ("unresolved facet orbit",d->pop!(d["representatives"][1]["neighbours"])),
        ("wrong transporter",d->begin
            changed=false
            for r in d["representatives"],n in r["neighbours"]
                if !changed && n["neighbour_gram"]!=d["representatives"][n["target"]]["gram"]
                    n["transporter_row"]=CONTRACT_AV.matrix_strings(identity_matrix(QQ,4));changed=true
                end
            end
            @assert changed
        end)]
    for (label,mutate) in mutations
        bad=deepcopy(good);payload=deepcopy(raw);mutate(payload)
        bad["raw_export_json"]=JSON3.write(payload)
        bad["raw_export_sha256"]=bytes2hex(CONTRACT_AV.sha256(bad["raw_export_json"]))
        file=joinpath(CONTRACT_ROOT,replace(label," "=>"_")*".json");CONTRACT_AV.atomic_json(file,bad)
        expect_backend_error(label,:arithmetic_validation) do;load_checkpoint(P,B,file);end
    end
    huge=ZZ(10)^80
    g=matrix(P.field,[1 P.field(huge)+P.generator;0 1]);R=encode_matrix(P,g)
    @test decode_matrix(P,R)==g
    @test CONTRACT_AV.matrix_wire(CONTRACT_AV.matrix_strings(R))==R
    OBSERVED["arbitrary precision wire roundtrip"]=Dict("status"=>"passed","decimal_digits"=>81)
end
CONTRACT_AV.atomic_json(joinpath(CONTRACT_ROOT,"summary.json"),
    Dict("status"=>"passed","diagnostic_root"=>CONTRACT_ROOT,"checks"=>OBSERVED))
println("Engineering and mathematical rejection contracts passed; diagnostics retained at ",CONTRACT_ROOT)
