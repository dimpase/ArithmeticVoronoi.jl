export export_dataset, verify_dataset, load_simplification

# Explicit archive I/O: verifies the recorded matrix maps, not a fresh
# abstract Tietze proof. Computation pipelines always simplify afresh.
function load_simplification(P,path,matrices,relations)
    d=JSON3.read(read(path,String),Dict{String,Any})
    S=(;label=d["group"],projective=get(d,"projective",false),algorithm=d["algorithm"],
      isomorphism_verified=d["isomorphism_verified"],old_generator_count=d["old_generator_count"],
      old_relator_count=d["old_relator_count"],gap_version=d["gap_version"],
      generator_matrices=Any[field_matrix_wire(P,g["matrix"]) for g in d["generators"]],
      backward_generator_words=[Int.(g["word_in_old_generators"]) for g in d["generators"]],
      forward_generator_words=[Int.(w) for w in d["old_generators_as_words_in_new"]],
      relators=[Int.(w) for w in d["relators"]])
    verify_simplification(matrices,relations,S)
    S
end



function export_dataset(results;directory=joinpath(pwd(),"results/dataset/v1"))
    mkpath(directory);artifacts=Dict{String,Any}()
    function copy_artifact(source,name)
        destination=joinpath(directory,name);mkpath(dirname(destination))
        abspath(source)==abspath(destination) || cp(source,destination;force=true)
        ref=Dict("path"=>name,"sha256"=>file_sha(destination),"bytes"=>filesize(destination))
        artifacts[name]=ref;ref
    end
    fields=Any[]
    for D in sort!(collect(keys(results)))
        R=results[D];P=R.enumeration.problem;prefix="fields/D$D"
        refs=Dict{String,Any}()
        for name in ["field.json","run.json","enumeration.json","GL_raw.json","SL_raw.json","PSL_raw.json",
                     "GL_simplified.json","SL_simplified.json","PSL_simplified.json","word_validation.json","initial_complete_run.json"]
            source=joinpath(R.directory,name)
            isfile(source) && (refs[name]=copy_artifact(source,"$prefix/$name"))
        end
        cache=joinpath(pwd(),"results/cache",R.enumeration.provenance["cache_key"]*".json")
        refs["checkpoint.json"]=copy_artifact(cache,"$prefix/checkpoint.json")
        sourcecert=normpath(joinpath(R.directory,R.run["certificate_directory"]))
        refs["coverage.json"]=copy_artifact(joinpath(sourcecert,"coverage.json"),"$prefix/two_cells/coverage.json")
        coverage=JSON3.read(read(joinpath(sourcecert,"coverage.json"),String),Dict{String,Any})
        for entry in coverage["entries"],file in entry["certificates"]
            copy_artifact(joinpath(sourcecert,file),"$prefix/two_cells/$file")
        end
        for sub in ["words","KW"]
            source=joinpath(R.directory,sub)
            isdir(source) || continue
            for file in sort!(readdir(source))
                endswith(file,".json") || continue
                copy_artifact(joinpath(source,file),"$prefix/$sub/$file")
            end
        end
        I=identity_matrix(P.field,2)
        minusGL=finite_element_word(R.gl,1,-I)
        embedding=[freely_reduce([repeat(R.sl.scalar_corrections[i]<0 ? inverse_word(R.sl.scalar_word) : R.sl.scalar_word,
            Int(abs(R.sl.scalar_corrections[i])));lift]) for (i,lift) in enumerate(R.sl.generator_lifts)]
        ensure(all(evaluate_word(R.gl.generator_matrices,w)==R.sl.generator_matrices[i] for (i,w) in enumerate(embedding)),"SL-to-GL word embedding failed")
        maps=(;schema="av_group_maps_v1",source_hash=R.context.source_hash,
            GL_scalar_quotient=(;domain="GL_raw",codomain="GL_mod_fundamental_scalar",
                generator_words=[[i] for i in eachindex(R.gl.generator_matrices)],additional_relator=R.sl.scalar_word,
                killed_scalar_coefficients=R.metadata.fundamental_unit,
                matrix_semantics="GL lifts modulo <eI>; central -I retained; this is not PGL"),
            SL_embedding_in_GL=(;domain="SL_raw",codomain="GL_raw",generator_words=embedding,
                method="prepend the scalar word to the recorded GL lift with its exact correction exponent"),
            determinant_square_kernel=(;domain="GL_mod_fundamental_scalar",quotient="C2 x C2",index=4,
                images=R.sl.determinant_images,transversal_GL_words=R.sl.transversal_words,
                schreier_ids=[collect(row) for row in eachrow(R.sl.schreier_ids)],
                meaning="squares in O^times; not a totally-positive-determinant condition"),
            SL_to_PSL=(;domain="SL_raw",codomain="PSL_raw",generator_words=[[i] for i in eachindex(R.sl.generator_matrices)],
                killed_central_word=R.sl.minus_identity_word,matrix_semantics="determinant-one lifts modulo +/-I"),
            minus_identity_GL_word=minusGL,minus_identity_SL_word=R.sl.minus_identity_word,
            fundamental_scalar_GL_word=R.sl.scalar_word)
        path=joinpath(R.directory,"group_maps.json");atomic_json(path,maps)
        refs["group_maps.json"]=copy_artifact(path,"$prefix/group_maps.json")
        groups=Any[]
        for name in ("GL","SL","PSL")
            arithmetic=(;field_discriminant=D,algebra="M2(K)",order="M2(O_K)",lattice="free O_K^2",
                rank=2,linearity="K-linear",semilinear_extension=false,
                determinant_condition=name=="GL" ? "unit" : "one",
                scalar_quotient=name=="PSL" ? "{I,-I}" : "none",
                contains_minus_identity=name!="PSL",includes_all_scalar_units=name=="GL")
            object_id="real_quadratic:D$D:Omax:free_rank2:K_linear:$name"
            raw=refs[name*"_raw.json"];simple=refs[name*"_simplified.json"]
            push!(groups,(;object_id,arithmetic,presentation_version="sha256:"*raw["sha256"],
                raw_presentation=raw,simplified_presentation=simple,
                matrix_semantics=name=="PSL" ? "SL lifts modulo +/-I" : "exact GL/SL matrices",
                invariants=R.run["presentations"][name]["abelianization"]))
        end
        push!(fields,(;discriminant=D,arithmetic_metadata=refs["field.json"],
            integral_basis_coordinates=[["1","0"],["0","1"]],groups,artifacts=refs,
            computation_run_id="sha256:"*refs["run.json"]["sha256"],
            geometry=R.run["geometry"],statuses=R.run["status"],
            certificate_directory="$prefix/two_cells",word_certificate_directory="$prefix/words",
            literature_directory=isfile(joinpath(R.directory,"KW/comparison.json")) ? "$prefix/KW" : nothing))
    end
    for file in readdir(joinpath(PACKAGE_ROOT,"src"))
        endswith(file,".jl") && copy_artifact(joinpath(PACKAGE_ROOT,"src",file),"software/src/$file")
    end
    for file in ["Project.toml","Manifest.toml"]
        copy_artifact(joinpath(PACKAGE_ROOT,file),"software/$file")
    end
    copy_artifact(joinpath(PACKAGE_ROOT,"scripts/verify_dataset.jl"),"software/scripts/verify_dataset.jl")
    for file in ["PRESENTATION_JUSTIFICATION.md","CONSTRUCTIVE_WORDS_STATUS.md","KW_COMPARISON.md","DATASET_SCHEMA.md","DATASET_PROPOSAL.md"]
        isfile(joinpath(PACKAGE_ROOT,file)) && copy_artifact(joinpath(PACKAGE_ROOT,file),"documentation/$file")
    end
    manifest=(;schema="arithmetic_groups_dataset_v1",dataset_version="local-v1",created_utc=string(now(UTC)),
        exact_encoding="field coefficients are rational strings in ascending powers of the recorded w; matrices are rows; words are signed one-based generator indices in written matrix-product order",
        provenance=(;workspace_git_commit=nothing,reason="workspace sources are not identified by a dedicated tracked commit; exact software hashes are supplied",
            backend_commit=let commits=unique([get(R.enumeration.provenance,"backend_commit",nothing) for R in values(results)])
                length(commits)==1 ? only(commits) : nothing
            end,
            backend_distribution="external pinned dependency; no upstream source trees or backend binaries bundled",
            julia_version=string(VERSION),oscar_version=string(pkgversion(Oscar)),
            hecke_version=string(pkgversion(Hecke)),json3_version=string(pkgversion(JSON3)),
            publication_status="local candidate dataset; not submitted or accepted by LMFDB; no DOI assigned"),
        fields,artifacts=artifacts,
        verification_levels=Dict("algebraic"=>"exact matrix membership, relators, maps, invariants, central elements and supplied path replay",
            "complete"=>"also independent raw perfect-enumeration validation, all ridge cycles, group-action assembly, SL/PSL extraction and fresh Tietze isomorphisms"))
    atomic_json(joinpath(directory,"dataset.json"),manifest)
    println("Dataset written to ",directory," (",length(artifacts)," hashed artifacts).");flush(stdout)
    directory
end

function check_projective_comparison(reference,computed,A,B,reference_relators,computed_relators)
    I=identity_matrix(base_ring(first(reference)),2)
    equivalent(x,y)=x==y || x==-y
    ensure(all(abs(k)<=length(reference) for w in B for k in w),"PSL comparison word uses a generator outside the reference alphabet")
    for (i,g) in enumerate(reference)
        ensure(equivalent(evaluate_word(computed,A[i]),g),"PSL reference-to-computed map failed")
        ensure(equivalent(evaluate_word(reference,substitute_word(A[i],B)),g),"PSL reference round trip failed")
    end
    for (i,g) in enumerate(computed)
        ensure(equivalent(evaluate_word(reference,B[i]),g),"PSL computed-to-reference map failed")
        ensure(equivalent(evaluate_word(computed,substitute_word(B[i],A)),g),"PSL computed round trip failed")
    end
    ensure(all(equivalent(evaluate_word(computed,substitute_word(w,A)),I) for w in reference_relators),"PSL substituted reference relation failed")
    ensure(all(equivalent(evaluate_word(reference,substitute_word(w,B)),I) for w in computed_relators),"PSL substituted computed relation failed")
    true
end

function dataset_path(root,reference)
    path=normpath(joinpath(root,reference isa AbstractString ? reference : reference["path"]))
    ensure(startswith(path,abspath(root)*"/"),"artifact path escapes the dataset directory")
    path
end
read_dataset_json(root,ref)=JSON3.read(read(dataset_path(root,ref),String),Dict{String,Any})

function wire_presentation(P,d)
    matrices=Any[field_matrix_wire(P,g["matrix"]) for g in d["generators"]]
    words=[Int.(w isa AbstractDict ? w["word"] : w) for w in d["relators"]]
    (;generator_matrices=matrices,relators=words)
end

function saved_word_context(ctx,d,objects,small,sd,units)
    P=ctx.enumeration.problem;graph=graph_edges(ctx);coherent=coherent_lifts(ctx,graph)
    ensure(field_matrix_data.(Ref(P),coherent.lifts)==d["coherent_vertex_lifts"],"saved coherent vertex lifts differ from the exact graph")
    matrices=objects[:GL].generator_matrices;groups=Any[];ids=Vector{Int}[]
    for g in d["vertex_groups"]
        localids=Int.(g["generator_ids"]);push!(ids,localids)
        elements=Any[field_matrix_wire(P,x) for x in g["elements"]]
        group=finite_group_record(P,elements;preferred=matrices[localids]);push!(groups,group)
        ensure(length(group.elements)==g["order"] && group.words==[Int.(w) for w in g["element_words"]],"saved finite stabilizer normal forms differ")
        ensure([collect(row) for row in eachrow(group.multiplication_table)]==g["multiplication_table"],"saved finite stabilizer multiplication table differs")
    end
    origins=Any[(;(Symbol(k)=>v for (k,v) in g["provenance"])...) for g in d["generators"]]
    C=PresentationCertificate(ctx,graph,coherent,groups,matrices,origins,objects[:GL].relators,Any[],ids,
        [Int.(h["word"]) for h in d["half_edges"]],
        [Int.(h["original_transporter_word"]) for h in d["half_edges"]],Int(d["attaching_loop_count"]))
    verify_presentation(C) # matrix consistency, not attaching-loop coverage
    idsdata=sd["schreier_ids"];schreier=[Int(idsdata[i][j]) for i in eachindex(idsdata),j in eachindex(first(idsdata))]
    S=SLPresentationCertificate(C,units,Int.(sd["scalar_GL_word"]),sd["scalar_word_search"],
        Int.(sd["GL_generator_determinant_mod_squares"]),[Int.(w) for w in sd["transversal_GL_words"]],schreier,
        [Int.(g["GL_lift_word"]) for g in sd["generators"]],
        [parse(BigInt,g["epsilon_scalar_exponent"]) for g in sd["generators"]],objects[:SL].generator_matrices,
        objects[:SL].relators,Any[],Int.(sd["minus_identity_word"]))
    verify_sl_presentation(S)
    C,S,psl_presentation(S)
end

function verify_dataset(directory;level=:algebraic,discriminants=nothing)
    root=abspath(directory);data=read_dataset_json(root,"dataset.json")
    ensure(data["schema"]=="arithmetic_groups_dataset_v1","unsupported dataset schema")
    level in (:algebraic,:complete) || throw(ArgumentError("level must be algebraic or complete"))
    for (name,ref) in data["artifacts"]
        ensure(name==ref["path"] && file_sha(dataset_path(root,ref))==ref["sha256"],"artifact checksum mismatch: $name")
    end
    observations=Any[]
    for field in data["fields"]
        D=Int(field["discriminant"])
        isnothing(discriminants) || D in discriminants || continue
        println("Verify D$D ($level)");flush(stdout)
        P=arithmetic_problem(D);I=identity_matrix(P.field,2);refs=field["artifacts"]
        metadata=read_dataset_json(root,refs["field.json"])
        ensure(metadata["polynomial_coefficients"]==[P.root_product,-P.root_sum,1] && discriminant(P.order)==D,"field/order basis mismatch")
        actual_metadata=field_metadata(P)
        ensure(actual_metadata.ordinary_class_number==metadata["ordinary_class_number"] && actual_metadata.fundamental_unit==metadata["fundamental_unit"],"class group or fundamental-unit data differs")
        objects=Dict{Symbol,Any}();small=Dict{Symbol,Any}();rawdata=Dict{Symbol,Any}()
        for name in (:GL,:SL,:PSL)
            d=read_dataset_json(root,refs[string(name)*"_raw.json"]);rawdata[name]=d
            obj=wire_presentation(P,d);objects[name]=obj
            for g in obj.generator_matrices
                check_group_matrix(P,encode_matrix(P,g));name==:GL || ensure(det(g)==1,"non-SL projective lift")
            end
            for (index,w) in enumerate(obj.relators)
                value=evaluate_word(obj.generator_matrices,w)
                ensure(value==I || (name==:PSL && value==-I),"raw $name relator failed")
                name==:PSL && ensure(value==Int(d["relators"][index]["SL_lift_sign"])*I,"PSL lift sign mismatch")
            end
            simple=load_simplification(P,dataset_path(root,refs[string(name)*"_simplified.json"]),obj.generator_matrices,obj.relators)
            small[name]=simple
            expected=only(g for g in field["groups"] if endswith(g["object_id"],":"*string(name)))["invariants"]
            actual=abelianization_invariants(length(simple.generator_matrices),simple.relators)
            ensure(actual.free_rank==expected["free_rank"] && actual.torsion==expected["torsion"],"abelianization mismatch")
        end
        units=verified_units(P);maps=read_dataset_json(root,refs["group_maps.json"])
        ensure(maps["GL_scalar_quotient"]["generator_words"]==[[i] for i in eachindex(objects[:GL].generator_matrices)] &&
            Int.(maps["GL_scalar_quotient"]["additional_relator"])==Int.(rawdata[:SL]["scalar_GL_word"]),"GL scalar-quotient map differs from its source alphabet")
        ensure(maps["SL_to_PSL"]["generator_words"]==[[i] for i in eachindex(objects[:SL].generator_matrices)] &&
            Int.(maps["SL_to_PSL"]["killed_central_word"])==Int.(rawdata[:SL]["minus_identity_word"]),"SL-to-PSL map differs from its source alphabet")
        ensure(evaluate_word(objects[:GL].generator_matrices,Int.(maps["minus_identity_GL_word"]))==-I,"GL lost -I")
        ensure(evaluate_word(objects[:SL].generator_matrices,Int.(maps["minus_identity_SL_word"]))==-I,"SL lost -I")
        ensure(evaluate_word(objects[:GL].generator_matrices,Int.(maps["fundamental_scalar_GL_word"]))==units.epsilon*I,"fundamental scalar word mismatch")
        for (i,w) in enumerate(maps["SL_embedding_in_GL"]["generator_words"])
            ensure(evaluate_word(objects[:GL].generator_matrices,Int.(w))==objects[:SL].generator_matrices[i],"SL/GL embedding map failed")
        end
        ensure(objects[:PSL].generator_matrices==objects[:SL].generator_matrices &&
            objects[:PSL].relators==[objects[:SL].relators;[Int.(maps["minus_identity_SL_word"])]] ,"PSL quotient map differs from its SL source")
        sd=rawdata[:SL];images=Int.(sd["GL_generator_determinant_mod_squares"])
        ensure(images==[unit_coordinates(P,units,det(g)).mod_squares for g in objects[:GL].generator_matrices],"determinant quotient images mismatch")
        ensure(length(sd["transversal_GL_words"])==4,"wrong Schreier quotient index")
        for (q,w) in enumerate(sd["transversal_GL_words"])
            ensure(foldl(xor,(images[abs(Int(j))] for j in w);init=0)==q-1,"Schreier transversal represents the wrong determinant coset")
        end
        for g in sd["generators"]
            lift=evaluate_word(objects[:GL].generator_matrices,Int.(g["GL_lift_word"]))
            ensure(field_matrix_wire(P,g["matrix"])==units.epsilon^parse(BigInt,g["epsilon_scalar_exponent"])*lift,"Schreier scalar correction mismatch")
        end
        # The existing validator consumes the exact exporter payload, without
        # invoking an executable or trusting serialized enumeration records.
        checkpoint=read_dataset_json(root,refs["checkpoint.json"])
        ensure(checkpoint["completed"] && checkpoint["problem"]==JSON3.read(JSON3.write(problem_descriptor(P)),Dict{String,Any}),"checkpoint field/action mismatch")
        ensure(bytes2hex(sha256(checkpoint["raw_export_json"]))==checkpoint["raw_export_sha256"],"checkpoint payload checksum mismatch")
        E=validate_export(P,JSON3.read(checkpoint["raw_export_json"],Dict{String,Any}))
        ctx=TwoCellContext(E)
        ensure(ctx.source_hash==rawdata[:GL]["source_hash"],"geometric source hash mismatch")
        coveragepath=joinpath(root,field["certificate_directory"],"coverage.json")
        coverage=JSON3.read(read(coveragepath,String),Dict{String,Any})
        cells=TwoCellCertificate[]
        if level==:complete
            for entry in coverage["entries"]
                pair=[read_two_cell(ctx,joinpath(dirname(coveragepath),file)) for file in entry["certificates"]]
                reverse_cycle_check(ctx,pair...);append!(cells,pair)
            end
            C=build_gl_presentation(ctx,cells;coverage)
            ensure(C.generator_matrices==objects[:GL].generator_matrices && C.relators==objects[:GL].relators,"reconstructed GL presentation differs from export")
            S=extract_sl_presentation(C);Q=psl_presentation(S)
            ensure(S.generator_matrices==objects[:SL].generator_matrices && S.relators==objects[:SL].relators,"reconstructed SL presentation differs from export")
            ensure(Q.relators==objects[:PSL].relators,"reconstructed PSL quotient differs from export")
            for name in (:GL,:SL,:PSL)
                fresh=simplify_matrix_presentation(objects[name].generator_matrices,objects[name].relators;label=small[name].label,projective=name==:PSL)
                ensure(fresh.forward_generator_words==small[name].forward_generator_words && fresh.backward_generator_words==small[name].backward_generator_words && fresh.relators==small[name].relators,"fresh Tietze replay differs; saved simplification is not re-certified")
            end
        else
            # Reconstruct the saved alphabets and coherent graph for path
            # replay, without accepting an attaching-loop coverage proof.
            C,S,Q=saved_word_context(ctx,rawdata[:GL],objects,small,sd,units)
        end
        R=ArithmeticGroupResult(E,ctx,C,S,Q,small,metadata,root,Dict{String,Any}())
        for path in readdir(joinpath(root,field["word_certificate_directory"]);join=true)
            endswith(path,".json") || continue
            cert=JSON3.read(read(path,String),Dict{String,Any});g=field_matrix_wire(P,cert["input_matrix"])
            verify_word_certificate(R,g,cert)
        end
        if !isnothing(field["literature_directory"])
            comp=read_dataset_json(root,field["literature_directory"]*"/comparison.json")
            refm=Any[field_matrix_wire(P,g) for g in comp["reference_matrices"]]
            A=[Int.(w) for w in comp["reference_to_computed_SL_words"]];B=[Int.(w) for w in comp["computed_SL_to_reference_words"]]
            ensure(all(evaluate_word(small[:SL].generator_matrices,w)==refm[i] for (i,w) in enumerate(A)),"literature forward map failed")
            ensure(all(evaluate_word(refm,w)==small[:SL].generator_matrices[i] for (i,w) in enumerate(B)),"literature backward map failed")
            PA=[Int.(w) for w in comp["reference_to_computed_PSL_words"]]
            PB=[Int.(w) for w in comp["computed_PSL_to_reference_words"]]
            check_projective_comparison(refm[1:4],small[:PSL].generator_matrices,PA,PB,
                [Int.(w) for w in comp["PSL_reference_relators"]],small[:PSL].relators)
            for (i,w) in enumerate(PB)
                ensure(evaluate_word(refm[1:4],w)==Int(comp["computed_PSL_to_reference_lift_signs"][i])*small[:PSL].generator_matrices[i],"PSL comparison lift sign differs")
                ensure(evaluate_word(refm,Int.(comp["computed_PSL_to_reference_SL_lift_words"][i]))==small[:PSL].generator_matrices[i],"PSL comparison exact SL-lift word differs")
            end
        end
        push!(observations,(;discriminant=D,status="passed",level=string(level)))
    end
    report=(;schema="av_dataset_verification_v1",status="passed",level=string(level),examples=observations,
        limitation=level==:algebraic ? "Algebraic consistency and path replay; no claim to re-certify complete two-cell coverage or saved Tietze isomorphisms at this level." : "Raw geometric completeness and fresh Tietze isomorphisms replayed; faithfulness invokes PRESENTATION_JUSTIFICATION.md, not matrix identities alone.")
    println(report.limitation);flush(stdout)
    report
end
