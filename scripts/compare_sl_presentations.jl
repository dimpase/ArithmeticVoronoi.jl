# Validation only: independent HAP/VMH presentations, exact maps, and Magma proofs.
include("compare_gl_presentations.jl")

function compare_sl_presentations(args;cases=((D=1,rank=3),(D=-4,rank=3),(D=1,rank=4),(D=-3,rank=3,abstract=false)))
    length(args) in (2,3) || error("provide real/imaginary VMH checkouts and optionally an output directory")
    roots=abspath.(args[1:2]);output=length(args)==3 ? abspath(args[3]) : abspath("results/sl_validation")
    mkpath(output)
    @assert GAP.Packages.load(get(ENV,"ARITHMETIC_VORONOI_HAP_ROOT","hap"))
    receipt=Dict{String,Any}("status"=>"running","julia_version"=>string(VERSION),
        "oscar_version"=>string(pkgversion(Oscar)),
        "hap_version"=>GAP.gap_to_julia(GAP.Globals.GAPInfo.PackagesLoaded.hap[2]),
        "vmh_commits"=>[readchomp(`git -C $root rev-parse HEAD`) for root in roots],
        "source_sha256"=>Dict(relpath(path,GPAV.PACKAGE_ROOT)=>bytes2hex(sha256(read(path)))
            for dir in ("src","gap") for path in readdir(joinpath(GPAV.PACKAGE_ROOT,dir);join=true) if isfile(path)),
        "cases"=>Any[])
    save()=GPAV.atomic_json(joinpath(output,"summary.json"),receipt)
    for case in cases
        case.D<=1 || error("fresh SL/PSL reference comparisons cover rational and imaginary quadratic fields")
        label="SL_PSL$(case.rank)_D$(case.D)";directory=joinpath(output,label);mkpath(directory)
        prove=get(case,:abstract,true)
        entry=Dict{String,Any}("rank"=>case.rank,"discriminant"=>case.D,"status"=>"running",
            "abstract_isomorphism_requested"=>prove)
        push!(receipt["cases"],entry);save()
        try
            P=case.D==1 ? classical_problem(rank=case.rank) : arithmetic_problem(case.D;rank=case.rank)
            R=run_field_pipeline(P;directory=joinpath(directory,"native"),progress=true)
            entry["native"]=R.run;entry["center_order"]=R.sl.center_order;save()
            println(label,": independent HAP SL presentation");flush(stdout)
            V=voronoi_complex(R.enumeration;group=:SL)
            F=GAP.Globals.FreeGResolution(hap_complex(V;resolution_length=2),2)
            T=GAP.Globals.PresentationOfResolution(F)
            data=GAP.gap_to_julia(GAP.evalstr("{R,P} -> [List(P.gens, i -> R!.elts[i]), List(P.relators, LetterRepAssocWord)]")(F,T);recursive=true)
            m=P.rank*P.degree;I=identity_matrix(P.field,P.rank)
            hm=[decode_matrix(P,matrix(QQ,m,m,[M[i][j] for i in 1:m for j in 1:m])) for M in data[1]]
            hr=[Int.(w) for w in data[2]]
            @assert all(w->GPAV.evaluate_word(hm,w)==I,hr)
            hc=[matrix_to_word(R,g;target=:SL) for g in hm]
            @assert all(verify_word_certificate(R,g,c) for (g,c) in zip(hm,hc))
            h1=sort(BigInt.(GAP.gap_to_julia(GAP.Globals.Homology(GAP.Globals.TensorWithIntegers(F),1))))
            @assert h1==presentation_abelianization(hm,hr)==presentation_abelianization(R.simplified[:SL].generator_matrices,R.simplified[:SL].relators)
            GPAV.atomic_json(joinpath(directory,"hap_sl.json"),(;matrices=GPAV.field_matrix_data.(Ref(P),hm),relators=hr,certificates=hc))
            entry["hap"]=Dict("SL_homology"=>h1,"SL_generators"=>length(hm),"SL_relators"=>length(hr),"matrix_relators_and_words"=>true);save()
            working=mktempdir(directory;prefix="vmh_",cleanup=false)
            root=roots[P.discriminant<0 ? 2 : 1]
            run(pipeline(`git -C $root archive HEAD`,`tar -x -C $working`))
            d=mod(case.D,4)==1 ? case.D : div(case.D,4)
            basic=case.D<0 ? "n:=$(case.rank); d:=$d; steinitz:=1;\n" :
                "K:=$(case.D==1 ? "QNF()" : "QuadraticField($d)"); n:=$(case.rank); V:=KMatrixSpace(K,n,1); L:=[IntegralBasis(K)[i]*V.j: i in [1..Degree(K)], j in [1..n]];\n"
            write(joinpath(working,"Magma/BasicData.m"),basic)
            write(joinpath(working,"parameters.m"),"FieldDiscriminant:=$(case.D); Dimension:=$(case.rank); Radicand:=$d;\n")
            cp(joinpath(@__DIR__,"magma_sl_presentations.m"),joinpath(working,"validate.m");force=true)
            println(label,": independent Magma SL and full-centre quotient");flush(stdout);start=time_ns()
            checked_magma(working,"validate.m",joinpath(directory,"magma.log"),"sl_complete.txt")
            reference=JSON3.read(read(joinpath(working,"reference_sl.json"),String),Dict{String,Any})
            @assert reference["center_order"]==R.sl.center_order
            entry["magma_seconds"]=(time_ns()-start)/1e9;entry["magma_version"]=reference["magma_version"]
            entry["magma"]=Dict{String,Any}();save()
            for target in (:SL,:PSL)
                ref=reference[string(target)];S=R.simplified[target]
                matrices=reference_matrix.(Ref(P),ref["matrices"])
                P.discriminant<0 && (matrices=[transpose(inv(g)) for g in matrices])
                relations=[Int.(w) for w in ref["relators"]]
                ntoR=short_generator_words(P,matrices,S.generator_matrices;projective=target==:PSL)
                equal(A,B)=A==B || (target==:PSL && !isnothing(GPAV.projective_scalar(A,B)))
                @assert all(g->det(g)==1 && check_group_matrix(P,encode_matrix(P,g)),matrices)
                @assert all(w->equal(GPAV.evaluate_word(matrices,w),I),relations)
                @assert all(equal(GPAV.evaluate_word(matrices,w),g) for (w,g) in zip(ntoR,S.generator_matrices))
                c=[matrix_to_word(R,g;target) for g in matrices]
                @assert all(verify_word_certificate(R,g,a) for (g,a) in zip(matrices,c))
                rtoN=short_generator_words(P,S.generator_matrices,matrices;fallback=[Int.(a["word"]) for a in c],radius=4,projective=target==:PSL)
                h1=presentation_abelianization(S.generator_matrices,S.relators)
                @assert h1==sort(BigInt.(ref["abelianization"]))==presentation_abelianization(matrices,relations)
                GPAV.atomic_json(joinpath(directory,"$(target)_maps.json"),(;reference_matrices=GPAV.field_matrix_data.(Ref(P),matrices),
                    reference_relators=relations,native_to_reference=ntoR,reference_to_native=rtoN,certificates=c))
                entry["magma"][string(target)]=Dict("abstract_isomorphism"=>(prove ? "pending" : "not_requested"),
                    "matrix_relators_and_bidirectional_words"=>true,"homology"=>h1,"generators"=>length(matrices),"relators"=>length(relations));save()
                if prove
                    write(joinpath(working,"presentation_maps.m"),join(["$name := $(JSON3.write(words));\n" for (name,words) in
                        (("NativeRelations",S.relators),("ReferenceRelations",relations),("NativeToReference",ntoR),("ReferenceToNative",rtoN))]))
                    cp(joinpath(@__DIR__,"magma_check_gl_maps.m"),joinpath(working,"maps.m");force=true)
                    println(label,": proving $target abstract isomorphism");flush(stdout)
                    checked_magma(working,"maps.m",joinpath(directory,"magma_$(target)_maps.log"),"maps_complete.txt")
                    entry["magma"][string(target)]["abstract_isomorphism"]="verified";save()
                end
            end
            entry["status"]="passed";save()
        catch err
            entry["status"]="failed";entry["error"]=sprint(showerror,err);receipt["status"]="failed";save();rethrow()
        end
    end
    receipt["status"]="passed";save();println("SL/PSL comparisons passed: ",joinpath(output,"summary.json"))
end
if abspath(PROGRAM_FILE)==@__FILE__
    compare_sl_presentations(ARGS)
end
