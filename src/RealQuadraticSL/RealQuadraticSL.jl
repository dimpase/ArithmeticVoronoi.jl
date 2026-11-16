module RealQuadraticSL

# The certified real-quadratic level construction used by the HAP provider.
include("morse_links.jl")
include("flag_links.jl")
include("level_complex.jl")
include("reduce_level.jl")

export compute_real_sl_complex

"""
    compute_real_sl_complex(D; rank=2, resolution_length=4, reduce=true, progress=false,
                            certificate_directory=nothing, max_cell_dimension=nothing,
                            hap_root="hap")

Compute a fresh, certified determinant-height level complex for SL(rank,O), where D
is a positive quadratic field discriminant. Returns `geometric`, its optional
algebraically reduced `chain`, and `hap` (a HAP non-free resolution), together with
link-collapse certificates and timings. `max_cell_dimension=m` explicitly requests
only the m-skeleton of the certified full level. For a proper skeleton,
`resolution_length` must not exceed m; homology is then supported through degree
`resolution_length-1`. The default constructs all dimensions.

Stabilizers retain the full scalar centre
of SL: -I is retained in even rank and excluded in odd rank.

The chosen rational height must pass exact ascending/descending
link collapse checks. Failure returns no complex; there is no precomputed-data
or homology-only fallback. Rank-two regression cases have discriminants 5 and 12.
The construction
provides no contracting homotopy for HAP cup-product computations.
"""
function compute_real_sl_complex(D::Integer;rank=2,kwargs...)
    D>1 || throw(ArgumentError("a positive real quadratic field discriminant is required"))
    compute_real_sl_complex(arithmetic_problem(Int(D);rank);kwargs...)
end

function compute_real_sl_complex(P::ArithmeticVoronoiProblem;resolution_length=4,reduce=true,
        progress=false,certificate_directory=nothing,max_cell_dimension=nothing,
        hap_root=get(ENV,"ARITHMETIC_VORONOI_HAP_ROOT","hap"))
    P.degree==2 && P.discriminant>1 || throw(ArgumentError("a real quadratic problem is required"))
    resolution_length isa Integer && resolution_length>0 || throw(ArgumentError("resolution_length must be positive"))
    if !isnothing(max_cell_dimension)
        max_cell_dimension isa Integer && !(max_cell_dimension isa Bool) && max_cell_dimension>=0 ||
            throw(ArgumentError("max_cell_dimension must be a nonnegative integer"))
        max_cell_dimension<P.rank^2-1 && resolution_length>max_cell_dimension &&
            throw(ArgumentError("a proper skeleton must include cells through resolution_length"))
    end
    D=P.discriminant
    timings=Dict{String,Float64}()
    function stage(label,f)
        if progress;println("Real SL",P.rank,": ",label);flush(stdout);end
        started=time_ns()
        value=f()
        timings[label]=(time_ns()-started)/1e9
        value
    end
    E=stage("fresh native enumeration",()->enumerate_perfect_forms(P;backend=:native,progress))
    C=stage("GL cell complex",()->voronoi_complex(E;progress))
    M=stage("cell orbit data",()->morse_experiment(C;progress,maximal_flags=false))
    data,links,faces=stage("ascending and descending links",()->morse_link_graphs(M;progress))
    certificate=stage("exact collapse certification",()->certify_morse_graphs(data;progress))
    L=stage("SL level cells",()->build_level_complex(M,data,certificate;progress,max_cell_dimension,faces))
    @assert validate_level_complex(L)
    reduced=reduce ? stage("equivariant chain cancellations",()->reduce_level_complex(L;progress)) : L
    @assert validate_level_complex(reduced)
    K=stage("HAP export",()->level_hap_complex(reduced;resolution_length,hap_root))
    if !isnothing(certificate_directory)
        mkpath(certificate_directory)
        write(joinpath(certificate_directory,"links.json"),JSON.json(data))
        witness=Dict("discriminant"=>D,"rank"=>P.rank,"link_model"=>"flag_graph","base_heights"=>string.(certificate.bases),
            "links"=>[Dict("orbit"=>r.orbit,"direction"=>r.sgn,
                "graph_collapse_steps"=>r.result.steps) for r in certificate.reports])
        write(joinpath(certificate_directory,"collapse_certificate.json"),JSON.json(witness))
        write(joinpath(certificate_directory,"construction.json"),JSON.json(Dict(
            "discriminant"=>D,"rank"=>P.rank,"actual_SL_complex"=>true,
            "certified_dimension"=>length(C.cells)-2,"computed_dimension"=>length(L.cells)-1,"cell_orbits"=>length.(L.cells),
            "chain_module_orbits"=>length.(reduced.cells),"timings"=>timings,
            "precomputed_data_used"=>false);pretty=4))
    end
    (;geometric=L,chain=reduced,hap=K,timings,certificate,M,data,links,L,reduced,K)
end
end
