module RealQuadraticSL

# Explicitly loaded prototype. The registered package/provider dispatch is unchanged.
include("morse_links.jl")
include("level_complex.jl")
include("reduce_level.jl")

export compute_real_sl_complex

"""
    compute_real_sl_complex(D; resolution_length=4, reduce=true, progress=false,
                            certificate_directory=nothing, hap_root="hap")

Compute a fresh, certified determinant-height level complex for SL₂(O), where D
is a positive quadratic field discriminant. Returns `geometric`, its optional
algebraically reduced `chain`, and `hap` (a HAP non-free resolution), together with
link-collapse certificates and timings. All stabilizers retain central -I.

Experimental: the chosen rational height must pass exact ascending/descending
link collapse checks. Failure returns no complex; there is no precomputed-data
or homology-only fallback. Tested on discriminants 5 and 12. The construction
provides no contracting homotopy for HAP cup-product computations.
"""
function compute_real_sl_complex(D::Integer;resolution_length=4,reduce=true,
        progress=false,certificate_directory=nothing,
        hap_root=get(ENV,"ARITHMETIC_VORONOI_HAP_ROOT","hap"))
    D>1 || throw(ArgumentError("a positive real quadratic field discriminant is required"))
    resolution_length isa Integer && resolution_length>0 || throw(ArgumentError("resolution_length must be positive"))
    P=arithmetic_problem(Int(D))
    timings=Dict{String,Float64}()
    function stage(label,f)
        if progress;println("Real SL₂: ",label);flush(stdout);end
        started=time_ns()
        value=f()
        timings[label]=(time_ns()-started)/1e9
        value
    end
    E=stage("fresh native enumeration",()->enumerate_perfect_forms(P;backend=:native,progress))
    C=stage("GL cell complex",()->voronoi_complex(E;progress))
    M=stage("barycentric flags",()->morse_experiment(C;progress))
    data,links=stage("ascending and descending links",()->morse_link_data(M))
    certificate=stage("exact collapse certification",()->certify_morse_links(data))
    L=stage("SL level cells",()->build_level_complex(M,data,certificate;progress))
    @assert validate_level_complex(L)
    reduced=reduce ? stage("equivariant chain cancellations",()->reduce_level_complex(L;progress)) : L
    @assert validate_level_complex(reduced)
    K=stage("HAP export",()->level_hap_complex(reduced;resolution_length,hap_root))
    if !isnothing(certificate_directory)
        mkpath(certificate_directory)
        write(joinpath(certificate_directory,"links.json"),JSON.json(data))
        witness=Dict("discriminant"=>D,"base_heights"=>string.(certificate.bases),
            "links"=>[Dict("orbit"=>r.orbit,"direction"=>r.sgn,
                "pairs"=>[[collect(a),collect(b)] for (a,b) in r.result.pairs]) for r in certificate.reports])
        write(joinpath(certificate_directory,"collapse_certificate.json"),JSON.json(witness))
        write(joinpath(certificate_directory,"construction.json"),JSON.json(Dict(
            "discriminant"=>D,"actual_SL_complex"=>true,"cell_orbits"=>length.(L.cells),
            "chain_module_orbits"=>length.(reduced.cells),"timings"=>timings,
            "precomputed_data_used"=>false);pretty=4))
    end
    (;geometric=L,chain=reduced,hap=K,timings,certificate,M,data,links,L,reduced,K)
end
end
