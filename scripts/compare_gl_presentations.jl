# Fresh, constructive GL comparisons; requires licensed Magma and both VMH checkouts.
using ArithmeticVoronoi, Oscar, JSON3, SHA
import Oscar: GAP
const GPAV=ArithmeticVoronoi

# For comparison, short generator maps make abstract relation proofs much
# smaller than the (valid but often very long) geometric descent words.
# Meet two balls in the exact integral matrix group; each returned word is
# checked. This is a validation utility, not the package's membership routine.
function short_generator_words(P,generators,targets;fallback=nothing,radius=typemax(Int))
    n=P.rank*P.degree
    matrices=encode_matrix.(Ref(P),generators)
    letters=[collect(eachindex(matrices));-collect(eachindex(matrices))]
    steps=[matrices;inv.(matrices)]
    values=[identity_matrix(QQ,n)];words=[Int[]]
    known=Dict(GPAV.matrix_key(first(values))=>1)
    encoded=encode_matrix.(Ref(P),targets)
    result=Union{Nothing,Vector{Int}}[nothing for _ in targets]
    frontier=1:1
    for depth in 0:radius
        for i in eachindex(targets)
            isnothing(result[i]) || continue
            for (j,g) in enumerate(values)
                k=get(known,GPAV.matrix_key(encoded[i]*inv(g)),0)
                k==0 && continue
                result[i]=GPAV.freely_reduce([words[k];words[j]])
                @assert GPAV.evaluate_word(generators,result[i])==targets[i]
                break
            end
        end
        println("  Generator-map search radius ",depth,": ",count(!isnothing,result),"/",length(result),
            " found; ",length(values)," matrices");flush(stdout)
        all(!isnothing,result) && return Vector{Int}[w for w in result]
        depth==radius && break
        firstnew=length(values)+1
        for j in frontier,(letter,g) in zip(letters,steps)
            value=values[j]*g;key=GPAV.matrix_key(value)
            haskey(known,key) && continue
            push!(values,value);push!(words,[words[j];letter]);known[key]=length(values)
        end
        frontier=firstnew:length(values)
        isempty(frontier) && error("reference generators close to a finite group without all requested matrices")
    end
    isnothing(fallback) && error("generator search did not find all requested matrices")
    [isnothing(w) ? Int.(fallback[i]) : w for (i,w) in enumerate(result)]
end

function reference_matrix(P,rows)
    root=P.degree==1 ? zero(P.field) : (P.root_sum==1 ? 2P.generator-1 : P.generator)
    # Entrywise field isomorphism: group matrices are not transposed.
    matrix(P.field,P.rank,P.rank,[sum(GPAV.rational_wire(c)*root^(k-1)
        for (k,c) in enumerate(rows[i][j])) for i in 1:P.rank for j in 1:P.rank])
end

function checked_magma(working,script,log,marker)
    rm(joinpath(working,marker);force=true)
    open(log,"w") do io
        run(pipeline(Cmd(`magma -b $script`;dir=working);stdout=io,stderr=io))
    end
    # Magma can return exit code zero after a script error.
    isfile(joinpath(working,marker)) || error("Magma did not finish: inspect $log")
end

function presentation_abelianization(matrices,relators)
    a=abelianization_invariants(length(matrices),relators)
    sort!([zeros(BigInt,a.free_rank);parse.(BigInt,a.torsion)])
end

function hap_presentation_comparison(R,directory)
    P=R.enumeration.problem;S=R.simplified[:GL]
    V=voronoi_complex(R.enumeration;progress=true)
    F=GAP.Globals.FreeGResolution(hap_complex(V;resolution_length=2),2)
    T=GAP.Globals.PresentationOfResolution(F)
    data=GAP.gap_to_julia(GAP.evalstr(
        "{R,P} -> [List(P.gens, i -> R!.elts[i]), List(P.relators, LetterRepAssocWord)]")(F,T);recursive=true)
    n=P.rank*P.degree
    matrices=[decode_matrix(P,matrix(QQ,n,n,[M[i][j] for i in 1:n for j in 1:n])) for M in data[1]]
    relators=[Int.(w) for w in data[2]]
    I=identity_matrix(P.field,P.rank)
    @assert all(w->GPAV.evaluate_word(matrices,w)==I,relators)
    certificates=[matrix_to_word(R,g) for g in matrices]
    @assert all(verify_word_certificate(R,g,c) for (g,c) in zip(matrices,certificates))
    words=[Int.(c["word"]) for c in certificates]
    h1=sort(BigInt.(GAP.gap_to_julia(GAP.Globals.Homology(GAP.Globals.TensorWithIntegers(F),1))))
    @assert h1==presentation_abelianization(S.generator_matrices,S.relators)==presentation_abelianization(matrices,relators)
    GPAV.atomic_json(joinpath(directory,"hap_presentation.json"),Dict(
        "matrices"=>GPAV.field_matrix_data.(Ref(P),matrices),"relators"=>relators,
        "words_in_native"=>words,"certificates"=>certificates))
    Dict("cells"=>length.(V.cells),"generators"=>length(matrices),"relators"=>length(relators),
        "abelianization"=>h1,"exact_relator_and_generator_word_checks"=>true)
end

function compare_gl_presentations(args;cases=((D=1,rank=3,hap=true),(D=-4,rank=3,hap=true),
        (D=1,rank=4,hap=true),(D=5,rank=3,hap=false),(D=-7,rank=2,hap=true)))
    length(args) in (2,3) || error("provide the real and imaginary VMH checkouts, and optionally an output directory")
    roots=abspath.(args[1:2])
    output=length(args)==3 ? abspath(args[3]) : joinpath(GPAV.PACKAGE_ROOT,"results/gl_validation")
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
        label="GL$(case.rank)_D$(case.D)"
        println("Comparing presentations: ",label);flush(stdout)
        directory=joinpath(output,label);mkpath(directory)
        entry=Dict{String,Any}("rank"=>case.rank,"discriminant"=>case.D,"status"=>"running")
        push!(receipt["cases"],entry);save()
        try
            P=case.D==1 ? classical_problem(;rank=case.rank) : arithmetic_problem(case.D;rank=case.rank)
            R=run_gl_pipeline(P;directory=joinpath(directory,"native"))
            S=R.simplified[:GL]
            entry["native"]=R.run;save()
            if case.hap
                println(label,": comparing HAP presentation and H1");flush(stdout)
                entry["hap"]=hap_presentation_comparison(R,directory);save()
            end
            # Copy every case separately. Preserve sources and old failed-run evidence.
            working=mktempdir(directory;prefix="vmh_",cleanup=false)
            rm(working);cp(roots[case.D<0 ? 2 : 1],working)
            d=mod(case.D,4)==1 ? case.D : div(case.D,4)
            basic=case.D<0 ? "n:=$(case.rank); d:=$d; steinitz:=1;\n" :
                "K:=$(case.D==1 ? "QNF()" : "QuadraticField($d)"); n:=$(case.rank); V:=KMatrixSpace(K,n,1); L:=[IntegralBasis(K)[i]*V.j: i in [1..Degree(K)], j in [1..n]];\n"
            write(joinpath(working,"Magma/BasicData.m"),basic)
            for name in ("reference.json","complete.txt","maps_complete.txt")
                rm(joinpath(working,name);force=true)
            end
            cp(joinpath(@__DIR__,"magma_gl_presentation.m"),joinpath(working,"validate.m");force=true)
            println(label,": fresh Magma presentation");flush(stdout);started=time_ns()
            checked_magma(working,"validate.m",joinpath(directory,"magma.log"),"complete.txt")
            entry["magma_seconds"]=(time_ns()-started)/1e9
            reference=JSON3.read(read(joinpath(working,"reference.json"),String))
            forms=[transpose(reference_matrix(P,f)) for f in reference.forms]
            forms=[f/shortest_data(trace_gram(P,f)).minimum for f in forms]
            # Imaginary VMH uses row matrices; real VMH already uses columns.
            # Apply g -> g^-t only for the imaginary implementation, then
            # align the exact perfect form and the native coherent vertex lift.
            matches=Int[]
            for f in forms
                candidates=findall(r->!isnothing(equivalence(P,f,r.form)),R.enumeration.representatives)
                @assert length(candidates)==1
                push!(matches,only(candidates))
            end
            @assert sort(matches)==collect(eachindex(R.enumeration.representatives))
            target=R.enumeration.representatives[first(matches)].form
            h=first(forms)==target ? identity_matrix(P.field,P.rank) : equivalence(P,first(forms),target)
            @assert !isnothing(h)
            # Native vertex groups use coherent lifts along the quotient tree,
            # not the bare stored form when its orbit is not the root vertex.
            h*=R.gl.coherent.lifts[first(matches)]
            images=reference_matrix.(Ref(P),reference.matrices)
            P.discriminant<0 && (images=[transpose(inv(g)) for g in images])
            matrices=[inv(h)*g*h for g in images]
            relators=[Int.(w) for w in reference.relators]
            @assert reference.perfect_forms==length(R.enumeration.representatives)
            @assert all(g->check_group_matrix(P,encode_matrix(P,g)),matrices)
            @assert all(w->GPAV.evaluate_word(matrices,w)==identity_matrix(P.field,P.rank),relators)
            h1=presentation_abelianization(S.generator_matrices,S.relators)
            @assert h1==sort(BigInt.(reference.abelianization))==presentation_abelianization(matrices,relators)
            ntoR=short_generator_words(P,matrices,S.generator_matrices)
            @assert length(ntoR)==length(S.generator_matrices)
            @assert all(GPAV.evaluate_word(matrices,w)==g for (w,g) in zip(ntoR,S.generator_matrices))
            certificates=[matrix_to_word(R,g) for g in matrices]
            @assert all(verify_word_certificate(R,g,c) for (g,c) in zip(matrices,certificates))
            rtoN=short_generator_words(P,S.generator_matrices,matrices;
                fallback=[Int.(c["word"]) for c in certificates],radius=5)
            GPAV.atomic_json(joinpath(directory,"generator_maps.json"),Dict(
                "reference_matrices"=>GPAV.field_matrix_data.(Ref(P),matrices),"reference_relators"=>relators,
                "native_to_reference"=>ntoR,"reference_to_native"=>rtoN,"certificates"=>certificates,
                "conjugating_matrix"=>GPAV.field_matrix_data(P,h),
                "reference_matrix_conversion"=>P.discriminant<0 ?
                    "h^-1 * transpose(g^-1) * h after the entrywise field isomorphism" :
                    "h^-1 * g * h after the entrywise field isomorphism",
                "map_method"=>"exact matrix meet-in-the-middle search; native geometric certificates also retained"))
            entry["magma"]=Dict("version"=>reference.magma_version,"perfect_forms"=>reference.perfect_forms,
                "generators"=>length(matrices),"relators"=>length(relators),"abelianization"=>h1,
                "exact_form_isometries"=>matches,
                "exact_relator_and_bidirectional_generator_word_checks"=>true,
                "abstract_isomorphism"=>"pending","working_directory"=>relpath(working,output))
            save()
            write(joinpath(working,"presentation_maps.m"),join(["$name := $(JSON3.write(words));\n"
                for (name,words) in (("NativeRelations",S.relators),("ReferenceRelations",relators),
                    ("NativeToReference",ntoR),("ReferenceToNative",rtoN))]))
            cp(joinpath(@__DIR__,"magma_check_gl_maps.m"),joinpath(working,"maps.m");force=true)
            println(label,": checking abstract presentation maps in Magma");flush(stdout);started=time_ns()
            checked_magma(working,"maps.m",joinpath(directory,"magma_maps.log"),"maps_complete.txt")
            entry["magma"]["abstract_isomorphism"]="verified"
            entry["magma"]["map_check_seconds"]=(time_ns()-started)/1e9
            entry["status"]="passed";save()
        catch err
            entry["status"]="failed";entry["error"]=sprint(showerror,err)
            receipt["status"]="failed";save();rethrow()
        end
    end
    receipt["status"]="passed";save()
    println("All GL presentation comparisons passed: ",joinpath(output,"summary.json"))
end

if abspath(PROGRAM_FILE)==@__FILE__
    compare_gl_presentations(ARGS)
end
