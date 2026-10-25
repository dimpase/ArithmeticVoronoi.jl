include(joinpath(@__DIR__,"../integrations/perfect_forms/adapter.jl"))
using .PerfectFormsAdapter
# Optional independent comparison with a freshly run C++ exporter.
# julia --project=integrations/perfect_forms test/native_comparison.jl [discriminants...]
using ArithmeticVoronoi, Oscar, Test
const NCAV=ArithmeticVoronoi

function compare_enumerations(native,reference)
    P=native.problem
    D=P.discriminant
    expected=Dict(1=>1,5=>2,8=>2,12=>3,13=>9,17=>34,21=>16,24=>22,28=>30,
        29=>39,33=>149,37=>97,40=>52,41=>144,44=>66,56=>144,60=>89)
    @test length(native.representatives)==length(reference.representatives)
    haskey(expected,D) && @test length(native.representatives)==expected[D]
    @test length(native.neighbours)==length(reference.neighbours)
    orbit_map=Int[]
    transporters=Any[]
    for r in native.representatives
        matches=Tuple{Int,Any}[]
        for (j,s) in enumerate(reference.representatives)
            det(r.trace_gram)==det(s.trace_gram) || continue
            g=equivalence(P,s.form,r.form)
            isnothing(g) || push!(matches,(j,g))
        end
        @test length(matches)==1
        j,g=only(matches)
        push!(orbit_map,j)
        push!(transporters,encode_matrix(P,g))
        s=reference.representatives[j]
        @test r.stabilizer_order==s.stabilizer_order
        @test length(r.facets)==length(s.facets)
        @test Set(NCAV.signed_key(NCAV.apply_integral(last(transporters),v)) for v in r.minimal_vectors)==Set(Tuple(v) for v in s.minimal_vectors)
    end
    @test sort(orbit_map)==collect(eachindex(reference.representatives))
    # Compare each crossed facet and its actual neighbouring form after changing
    # source coordinates and reducing by the reference stabilizer. Counts alone
    # would not detect incorrect transporter or facet-orbit labels.
    for edge in native.neighbours
        source=reference.representatives[orbit_map[edge.source]]
        index=Dict(Tuple(v)=>k for (k,v) in enumerate(source.minimal_vectors))
        R=transporters[edge.source]
        local_vectors=native.representatives[edge.source].minimal_vectors
        mapped=sort([index[NCAV.signed_key(NCAV.apply_integral(R,local_vectors[k]))] for k in edge.facet])
        facet=only(f for f in source.facets if f.indices==mapped)
        edges=[e for e in reference.neighbours if e.source==orbit_map[edge.source]]
        other=edges[facet.neighbour_orbit]
        @test other.target==orbit_map[edge.target]
        S=encode_matrix(P,facet.to_orbit_representative)
        @test transpose(inv(R))*edge.neighbour_gram*inv(R)==transpose(S)*other.neighbour_gram*S
    end
    orbit_map
end

function compare_native(D,directory)
    P=D==1 ? classical_problem() : arithmetic_problem(D)
    native=enumerate_perfect_forms(P;backend=:native,progress=true)
    reference=enumerate_perfect_forms(P;backend=:polyhedral_common,force=true,
        cache_dir=joinpath(directory,"cache"),work_dir=joinpath(directory,"runs"),timeout=3600.0)
    compare_enumerations(native,reference)
    println("D$D comparison passed: ",length(native.representatives)," orbits, ",length(native.neighbours),
        " labelled neighbours; native seconds ",native.provenance["backend_seconds"],
        ", C++ seconds ",reference.provenance["backend_seconds"],
        ", C++ import/validation seconds ",reference.provenance["validation_seconds"])
end

if abspath(PROGRAM_FILE)==@__FILE__
    @testset "native/C++ exact comparison" begin
        mktempdir() do directory
            for D in (isempty(ARGS) ? [1,5,8,12,13,17,21,24,28,29,33,37,40,41,44,56,60] : parse.(Int,ARGS))
                @testset "D$D" begin
                    compare_native(D,directory)
                end
            end
        end
    end
end
