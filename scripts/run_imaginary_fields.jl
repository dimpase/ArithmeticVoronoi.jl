# From the checkout root: julia --project=. scripts/run_imaginary_fields.jl
# Supports --output=DIR, --resume, and explicit negative field discriminants.
include("run_native_fields.jl")
using SHA

const IMAGINARY_DISCRIMINANTS=[-3,-4,-7,-15,-23]
# Yasaki, arXiv:0908.1762v2, Tables 2--4; tuples count vertices, edges, faces.
const YASAKI_POLYTOPES=Dict(-3=>[(4,6,4)],-4=>[(6,12,8)],-7=>[(6,9,5)],
    -15=>[(4,6,4),(6,12,8)],-23=>[(5,8,5),(6,9,5),(6,12,8)])
const PAGE_BENCHMARKS=Dict(-3=>(0.93,0.83),-15=>(0.92,2.32),-23=>(1.22,2.06))

function imaginary_checks(R)
    P=R.enumeration.problem;D=P.discriminant
    D<0 || error("the imaginary survey requires negative field discriminants")
    shapes=sort!([(length(unique(NAV.primitive_ray(v) for v in NAV.matrix_rows(r.evaluation_rays))),
        length(R.context.ridges[i]),length(r.facets)) for (i,r) in enumerate(R.enumeration.representatives)])
    if haskey(YASAKI_POLYTOPES,D)
        NAV.ensure(shapes==YASAKI_POLYTOPES[D],"Voronoi polytope types disagree with Yasaki")
    end
    # Independent HAP complexes, distributed with Oscar's GAP installation.
    # Its input is a squarefree radicand, whereas our API uses discriminants.
    d=mod(D,4)==0 ? div(D,4) : D
    hap=Dict{String,Any}()
    for name in (:SL,:PSL)
        resolution=(name==:SL ? Oscar.GAP.Globals.ResolutionSL2QuadraticIntegers :
            Oscar.GAP.Globals.ResolutionPSL2QuadraticIntegers)(d,2)
        resolution==Oscar.GAP.Globals.fail && error("HAP has no comparison complex for $D")
        invariants=Int.(Oscar.GAP.gap_to_julia(Oscar.GAP.Globals.Homology(
            Oscar.GAP.Globals.TensorWithIntegers(resolution),1)))
        expected=R.run["presentations"][string(name)]["abelianization"]
        NAV.ensure(count(iszero,invariants)==expected.free_rank &&
            string.(filter(x->x>1,invariants))==expected.torsion,"$name abelianization disagrees with HAP")
        hap[string(name)]=Dict("status"=>"matched","invariant_factors"=>invariants)
    end
    # Timings use fresh enumerations: one warm-up, then three measured runs.
    enumerate_perfect_forms(P;backend=:native)
    times=Float64[]
    for _ in 1:3
        GC.gc();start=time_ns();E=enumerate_perfect_forms(P;backend=:native)
        push!(times,(time_ns()-start)/1e9)
        NAV.ensure(length(E.representatives)==length(R.enumeration.representatives),"fresh enumeration orbit count changed")
    end
    comparison=Dict{String,Any}("status"=>"passed",
        "polytope_f_vectors"=>collect.(shapes),
        "yasaki"=>Dict("status"=>haskey(YASAKI_POLYTOPES,D) ? "matched" : "not_tabulated_here",
            "source"=>"https://arxiv.org/abs/0908.1762","tables"=>"2, 3, 4"),
        "HAP"=>Dict("version"=>String(Oscar.GAP.gap_to_julia(Oscar.GAP.Globals.PackageInfo(Oscar.GAP.GapObj("hap"))[1].Version)),
            "source"=>"https://gap-packages.github.io/hap/tutorial/chap14.html",
            "radicand"=>d,"groups"=>hap,
            "scope"=>"independent SL/PSL abelianization comparison; not a certified isomorphism of presentations"),
        "native_enumeration_seconds"=>times,"native_enumeration_median_seconds"=>sort(times)[2])
    if haskey(PAGE_BENCHMARKS,D)
        psl,pgl=PAGE_BENCHMARKS[D]
        comparison["Page"]=Dict("source"=>"https://www.normalesup.org/~page/Recherche/Documents/articles/kln_gps.pdf",
            "table"=>"3.1","published_KleinianGroups_PSL_seconds"=>psl,
            "published_KleinianGroups_PGL_seconds"=>pgl,
            "status"=>"historical published timings; Magma was not rerun",
            "scope"=>"different hardware and tasks; these timings do not establish a speed ratio")
    end
    NAV.atomic_json(joinpath(R.directory,"comparison.json"),comparison)
    R.run["status"]["literature_comparison"]="Yasaki_polytope_types_and_HAP_abelianizations_checked"
    NAV.atomic_json(joinpath(R.directory,"run.json"),R.run)
    comparison
end

function imaginary_bundle(root,summary)
    readjson(path)=JSON3.read(read(path,String),Dict{String,Any})
    fields=Dict{String,Any}()
    for D in summary["discriminants"]
        directory=joinpath(root,"D$D");entry=summary["fields"][string(D)]
        groups=Dict{String,Any}()
        for name in ("GL","SL","PSL")
            path=joinpath(directory,"$(name)_simplified.json");saved=readjson(path)
            groups[name]=Dict("generator_matrices"=>[g["matrix"] for g in saved["generators"]],
                "relators"=>saved["relators"],"projective"=>name=="PSL",
                "abelianization"=>entry["presentations"][name]["abelianization"],
                "certificate_sha256"=>bytes2hex(sha256(read(path))))
        end
        fields[string(D)]=Dict("field"=>readjson(joinpath(directory,"field.json")),
            "perfect_form_orbits"=>entry["perfect_form_orbits"],
            "labelled_neighbours"=>entry["labelled_neighbours"],"groups"=>groups,
            "comparison"=>entry["comparison"],"timings"=>entry["stage_seconds"],
            "pipeline_seconds"=>entry["pipeline_seconds"],"word_validation"=>entry["word_validation"])
    end
    bundle=Dict("schema"=>"av_imaginary_presentations_v1","status"=>"complete",
        "discriminants"=>summary["discriminants"],"julia_version"=>string(VERSION),
        "oscar_version"=>string(pkgversion(Oscar)),"julia_threads"=>Threads.nthreads(),
        "blas_threads"=>summary["blas_threads"],"timing_note"=>summary["timing_note"],
        "enumeration_timing_note"=>"One warm-up followed by three fresh native enumerations per field, with GC before each; package import and compilation excluded.",
        "cpu"=>Sys.cpu_info()[1].model,
        "matrix_encoding"=>"Each entry [a,b] means a+b*w, with exact rational strings and w=(s+sqrt(D))/2, s=D mod 4.",
        "word_encoding"=>"Signed one-based generator indices in matrix-product order. PSL matrices are SL lifts and relators evaluate to +/-I.",
        "certificates"=>"Full raw presentations, Tietze generator maps, ridge coverage and word certificates are in the run directory; this bundle stores the simplified presentations.",
        "fields"=>fields)
    NAV.atomic_json(joinpath(root,"presentations.json"),bundle)
    bundle
end

if abspath(PROGRAM_FILE)==@__FILE__
    Oscar.GAP.Packages.load("hap") || error("the comparison script requires GAP's HAP package")
    result=main(ARGS;on_result=imaginary_checks,default_discriminants=IMAGINARY_DISCRIMINANTS,
        default_output=joinpath(NAV.PACKAGE_ROOT,"results/imaginary_fields"))
    imaginary_bundle(result.root,result.summary)
end
