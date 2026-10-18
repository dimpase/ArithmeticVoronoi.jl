using ArithmeticVoronoi, Oscar, JSON3, SHA, Dates
const AV=ArithmeticVoronoi
const ROOT=normpath(joinpath(@__DIR__,".."))

function word_text(word;prefix="g")
    isempty(word) && return "1"
    terms=String[];i=1
    while i<=length(word)
        j=i+1
        while j<=length(word) && word[j]==word[i];j+=1;end
        exponent=sign(word[i])*(j-i)
        push!(terms,prefix*string(abs(word[i]))*(exponent==1 ? "" : "^"*string(exponent)))
        i=j
    end
    join(terms," ")
end

function report_one(io,result)
    println(io,"Discriminant ",result["discriminant"],": ",result["status"])
    for key in ["local_ridge_counts","locally_reduced_ridge_counts","oriented_loops","loop_length_counts",
                "nonidentity_residuals","minus_identity_residuals","global_unoriented_edges","inverted_edges",
                "GL_raw_generators","GL_raw_relators","SL_raw_generators","SL_raw_relators",
                "GL_simplified_generators","GL_simplified_relators","SL_simplified_generators","SL_simplified_relators",
                "stage_seconds"]
        haskey(result,key) && println(io,"  ",key,": ",result[key])
    end
end

function run_case(D)
    directory=joinpath(ROOT,"results/two_cells","D$D");mkpath(directory)
    times=Dict{String,Float64}()
    times["validated_enumeration_load"]=@elapsed E=enumerate_perfect_forms(arithmetic_problem(D))
    times["ridge_context"]=@elapsed ctx=TwoCellContext(E)
    times["loops_and_coverage"]=@elapsed cells=enumerate_two_cells(ctx;directory)
    write_first_two_cell(ctx,first(cells.certificates);directory)
    times["GL_assembly_and_validation"]=@elapsed C=build_gl_presentation(ctx,cells.certificates;coverage=cells.ledger)
    atomic=AV.atomic_json
    raw=JSON3.read(JSON3.write(presentation_data(C)),Dict{String,Any})
    raw["two_cell_certificates"]=[Dict("source"=>cell.source,"ridge"=>cell.ridge,"initial_outgoing"=>cell.initial_outgoing,
        "sha256"=>bytes2hex(sha256(JSON3.write(two_cell_data(ctx,cell))))) for cell in cells.certificates]
    raw["coverage_file_sha256"]=AV.file_sha(joinpath(directory,"coverage.json"))
    atomic(joinpath(directory,"GL_presentation.json"),raw)
    times["SL_extraction_and_validation"]=@elapsed S=extract_sl_presentation(C)
    atomic(joinpath(directory,"SL_presentation.json"),sl_presentation_data(S))
    times["GL_Tietze_and_maps"]=@elapsed Gsmall=simplify_matrix_presentation(C.generator_matrices,C.relators;label="GL(2,O), D$D")
    atomic(joinpath(directory,"GL_simplified.json"),AV.simplification_data(E.problem,Gsmall))
    times["SL_Tietze_and_maps"]=@elapsed Ssmall=simplify_matrix_presentation(S.generator_matrices,S.relators;label="SL(2,O), D$D")
    atomic(joinpath(directory,"SL_simplified.json"),AV.simplification_data(E.problem,Ssmall))
    identity=identity_matrix(E.problem.field,2)
    first_nonidentity=findfirst(cell->cell.residual!=identity,cells.certificates)
    if !isnothing(first_nonidentity)
        atomic(joinpath(directory,"nontrivial_residual_two_cell.json"),two_cell_data(ctx,cells.certificates[first_nonidentity]))
    end
    lengths=Dict{String,Int}()
    for cell in cells.certificates
        key=string(length(cell.steps));lengths[key]=get(lengths,key,0)+1
    end
    result=Dict{String,Any}("discriminant"=>D,"status"=>"certified_GL_and_SL","source_hash"=>ctx.source_hash,
        "local_ridge_counts"=>length.(ctx.ridges),"locally_reduced_ridge_counts"=>cells.ledger["locally_reduced_ridge_counts"],
        "oriented_loops"=>length(cells.certificates),"loop_length_counts"=>lengths,
        "nonidentity_residuals"=>count(cell->cell.residual!=identity,cells.certificates),
        "minus_identity_residuals"=>count(cell->cell.residual==-identity,cells.certificates),
        "global_unoriented_edges"=>length(C.graph.edges),"inverted_edges"=>count(e->!isempty(e.reversals),C.graph.edges),
        "GL_raw_generators"=>length(C.generator_matrices),"GL_raw_relators"=>length(C.relators),
        "SL_raw_generators"=>length(S.generator_matrices),"SL_raw_relators"=>length(S.relators),
        "GL_simplified_generators"=>length(Gsmall.generator_matrices),"GL_simplified_relators"=>length(Gsmall.relators),
        "SL_simplified_generators"=>length(Ssmall.generator_matrices),"SL_simplified_relators"=>length(Ssmall.relators),
        "scalar_unit_word"=>S.scalar_word,"scalar_word_search"=>S.scalar_reduction,
        "stage_seconds"=>times,"loop_cache_hits"=>cells.ledger["cache_hits"],
        "source_enumeration_cache_key"=>E.provenance["cache_key"],
        "global_two_cell_orbits"=>nothing,"global_two_cell_reduction"=>"not performed; redundant exhaustive local coverage",
        "proof"=>"PRESENTATION_JUSTIFICATION.md","presentation_directory"=>relpath(directory,ROOT))
    atomic(joinpath(directory,"summary.json"),result)
    open(joinpath(directory,"PRESENTATIONS.txt"),"w") do io
        println(io,"Exact arithmetic Voronoi presentations; D=",D)
        println(io,"GL(2,O) and SL(2,O), K-linear; scalar units included in GL, -I retained in both.")
        println(io,"w satisfies ",D==5 ? "w^2-w-1=0" : "w^2-2=0","; all matrices use this maximal-order basis.")
        report_one(io,result)
        for (label,presentation) in [("GL",Gsmall),("SL",Ssmall)]
            println(io,"\n",label," simplified generators:")
            for (i,g) in enumerate(presentation.generator_matrices);println(io,"  g",i," = ",g);end
            println(io,label," relators (written matrix-product order):")
            for (i,word) in enumerate(presentation.relators);println(io,"  r",i,": ",word_text(word)," = 1");end
        end
        minus_gl=AV.finite_element_word(C,1,-identity)
        println(io,"\n-I in simplified GL generators: ",word_text(AV.substitute_word(minus_gl,Gsmall.forward_generator_words)))
        println(io,"-I in simplified SL generators: ",word_text(AV.substitute_word(S.minus_identity_word,Ssmall.forward_generator_words)))
        println(io,"Both Tietze maps are saved in both directions in *_simplified.json.")
        println(io,"The central scalar quotient and determinant-one rescaling are saved in SL_presentation.json.")
        println(io,"Completeness justification: ../../../PRESENTATION_JUSTIFICATION.md (at the package root).")
    end
    report_one(stdout,result);flush(stdout)
    result
end

selected=isempty(ARGS) ? [5,8] : parse.(Int,ARGS)
report=Dict{String,Any}("schema"=>"arithmetic_voronoi_two_cell_summary_v1","julia_version"=>string(VERSION),
    "oscar_version"=>string(pkgversion(Oscar)),"started_utc"=>string(now(UTC)),"examples"=>Any[])
for D in selected
    push!(report["examples"],run_case(D))
    AV.atomic_json(joinpath(ROOT,"results/two_cell_summary.json"),report)
end
