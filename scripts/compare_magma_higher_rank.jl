# julia --project=. scripts/compare_magma_higher_rank.jl /path/to/VMH-TotallyRealNumberFields /path/to/VMH-ImaginaryQuadraticNumberFields [output-directory] [--extended]
# Requires licensed Magma. Each case copies VMH and changes only BasicData.m;
# the supplied source checkouts and all upstream algorithms are left intact.
using ArithmeticVoronoi, Oscar, JSON3, SHA
import Oscar: GAP

function compare_magma_higher_rank(args;cases=((D=5,rank=2,complex=true),
        (D=-7,rank=2,complex=true),(D=1,rank=3,complex=true),
        (D=1,rank=4,complex=false),(D=5,rank=3,complex=false),
        (D=-4,rank=3,complex=true,homology_degree=1)))
    extended="--extended" in args
    args=filter(!=("--extended"),args)
    length(args) in (2,3) || error("provide the real and imaginary VMH checkouts, and optionally an output directory")
    roots=abspath.(args[1:2])
    output=length(args)==3 ? abspath(args[3]) : joinpath(ArithmeticVoronoi.PACKAGE_ROOT,"results/higher_rank_validation")
    mkpath(output)
    hap_root=get(ENV,"ARITHMETIC_VORONOI_HAP_ROOT","hap")
    @assert GAP.Packages.load(hap_root)
    function homology(R,degree)
        C=GAP.Globals.TensorWithIntegers(R)
        println("  Integral chain ranks: ",GAP.gap_to_julia(GAP.evalstr(
            "{C,n} -> List([0..n], C!.dimension)")(C,degree+1)));flush(stdout)
        result=Vector{BigInt}[]
        for n in 1:degree
            println("  Computing integral H",n);flush(stdout)
            started=time_ns()
            push!(result,sort(BigInt.(GAP.gap_to_julia(GAP.Globals.Homology(C,n)))))
            println("  H",n," = ",last(result),"; seconds=",(time_ns()-started)/1e9);flush(stdout)
        end
        result
    end
    receipt=Dict{String,Any}("status"=>"running","julia_version"=>string(VERSION),
        "oscar_version"=>string(pkgversion(Oscar)),
        "hap_version"=>GAP.gap_to_julia(GAP.Globals.GAPInfo.PackagesLoaded.hap[2]),
        "vmh_commits"=>[readchomp(`git -C $root rev-parse HEAD`) for root in roots],
        "provider_source_sha256"=>Dict(relpath(path,ArithmeticVoronoi.PACKAGE_ROOT)=>bytes2hex(sha256(read(path)))
            for dir in ("src","gap") for path in readdir(joinpath(ArithmeticVoronoi.PACKAGE_ROOT,dir);join=true)
            if isfile(path)),
        "cases"=>Any[])
    save()=ArithmeticVoronoi.atomic_json(joinpath(output,"summary.json"),receipt)
    for case in cases
        label="rank$(case.rank)_D$(case.D)"
        println("Comparing ",label);flush(stdout)
        degree=extended ? 2 : get(case,:homology_degree,2)
        entry=Dict{String,Any}("rank"=>case.rank,"discriminant"=>case.D,"status"=>"running",
            "homology_degree"=>case.complex ? degree : 0)
        push!(receipt["cases"],entry);save()
        try
            root=roots[case.D<0 ? 2 : 1]
            directory=mktempdir(output;prefix=label*"_",cleanup=false)
            working=joinpath(directory,"vmh")
            cp(root,working)
            d=mod(case.D,4)==1 ? case.D : div(case.D,4)
            basic=if case.D<0
                "n:=$(case.rank); d:=$d; steinitz:=1;\n"
            else
                field=case.D==1 ? "QNF()" : "QuadraticField($d)"
                "K:=$field; n:=$(case.rank); V:=KMatrixSpace(K,n,1); L:=[IntegralBasis(K)[i]*V.j: i in [1..Degree(K)], j in [1..n]];\n"
            end
            write(joinpath(working,"Magma/BasicData.m"),basic)
            # VMH's Write appends: copied files from an earlier local run must
            # never become part of the newly computed comparison.
            for name in ("forms.json","abelianization.json","complex_gl.g","complex_sl.g","complete.txt","magma_version.json")
                rm(joinpath(working,name);force=true)
            end
            write(joinpath(working,"validate.m"),"AVComputeComplex := $(case.complex);\n"*
                read(joinpath(@__DIR__,"magma_voronoi_comparison.m"),String))
            started=time_ns()
            open(joinpath(directory,"magma.log"),"w") do io
                run(pipeline(Cmd(`magma -b validate.m`;dir=working);stdout=io,stderr=io))
            end
            @assert isfile(joinpath(working,"complete.txt")) "Magma did not finish; inspect $(joinpath(directory,"magma.log"))"
            entry["magma_seconds"]=(time_ns()-started)/1e9
            entry["magma_version"]=JSON3.read(read(joinpath(working,"magma_version.json"),String))
            entry["directory"]=relpath(directory,output)
            P=case.D==1 ? classical_problem(;rank=case.rank) : arithmetic_problem(case.D;rank=case.rank)
            started=time_ns()
            E=enumerate_perfect_forms(P;backend=:native,progress=true)
            entry["enumeration_seconds"]=(time_ns()-started)/1e9
            forms=JSON3.read(read(joinpath(working,"forms.json"),String))
            @assert length(E.representatives)==length(forms)
            root_image=P.degree==1 ? zero(P.field) : (P.root_sum==1 ? 2P.generator-1 : P.generator)
            matches=Int[]
            for rows in forms
                F=matrix(P.field,P.rank,P.rank,[sum(ArithmeticVoronoi.rational_wire(c)*root_image^(k-1)
                    for (k,c) in enumerate(rows[i][j])) for i in 1:P.rank for j in 1:P.rank])
                # VMH uses row vectors and x F conjugate(x)^t.
                F=transpose(F)
                F/=shortest_data(trace_gram(P,F)).minimum
                candidates=findall(r->!isnothing(equivalence(P,r.form,F)),E.representatives)
                @assert length(candidates)==1
                push!(matches,only(candidates))
            end
            @assert sort(matches)==collect(eachindex(E.representatives))
            entry["perfect_forms"]=length(forms)
            entry["exact_form_isometries"]=matches
            entry["stabilizer_orders"]=[r.stabilizer_order for r in E.representatives]
            entry["groups"]=Dict{String,Any}()
            save()
            if case.complex
                GAP.Globals.Read(GAP.GapObj(joinpath(working,"GAP/ReadWellRoundedComplex.gi")))
                for group in (case.D<0 ? (:GL,:SL) : (:GL,))
                    started=time_ns()
                    C=voronoi_complex(E;group,progress=true)
                    println("Native ",group," cells: ",length.(C.cells));flush(stdout)
                    native=homology(GAP.Globals.FreeGResolution(hap_complex(C;resolution_length=degree+1),degree+1),degree)
                    path=joinpath(working,"complex_$(lowercase(string(group))).g")
                    K=GAP.Globals.ResolutionFromWellRoundedComplex(GAP.GapObj(path),degree+1)
                    dims=Int.(GAP.gap_to_julia(GAP.Globals.DIMS))
                    @assert dims==length.(C.cells)
                    for k in eachindex(C.cells)
                        sizes=Int.(GAP.gap_to_julia(GAP.evalstr(
                            "{K,d,n} -> List([1..n], j -> Size(K!.stabilizer(d,j)))")(K,k-1,dims[k])))
                        @assert sort(sizes)==sort(length.(C.stabilizers[k]))
                    end
                    println("Magma ",group," comparison");flush(stdout)
                    reference=homology(GAP.Globals.FreeGResolution(K,degree+1),degree)
                    @assert native==reference
                    if group==:GL
                        abelianization=sort(BigInt.(JSON3.read(read(joinpath(working,"abelianization.json"),String))))
                        @assert native[1]==abelianization
                    end
                    result=Dict("cells"=>dims,"homology"=>native,"seconds"=>(time_ns()-started)/1e9,
                        "magma_complex_sha256"=>bytes2hex(sha256(read(path))))
                    if group==:GL && case.rank==3 && case.D<0
                        stored=homology(GAP.Globals.ResolutionGL3QuadraticIntegers(d,degree+1),degree)
                        @assert native==stored
                        result["hap_stored_homology"]=stored
                    end
                    entry["groups"][string(group)]=result;save()
                end
            end
            entry["status"]="passed";save()
        catch err
            entry["status"]="failed";entry["error"]=sprint(showerror,err)
            receipt["status"]="failed";save();rethrow()
        end
    end
    receipt["status"]="passed";save()
    println("All comparisons passed: ",joinpath(output,"summary.json"))
end

if abspath(PROGRAM_FILE)==@__FILE__
    compare_magma_higher_rank(ARGS)
end
