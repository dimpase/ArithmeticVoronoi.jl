using ArithmeticVoronoi, JSON3, SHA, Test
using .PerfectFormsAdapter
const COAV=ArithmeticVoronoi

@testset "C++ checkpoint loading requires explicit input" begin
    mktempdir() do directory
        # An exporter that deliberately fails lets us distinguish an actual
        # launch from replay, without requiring an installed C++ toolchain.
        executable=joinpath(directory,"exporter")
        write(executable,"#!/bin/sh\nexit 37\n");chmod(executable,0o700)
        B=BackendConfig(executable);P=arithmetic_problem(5)
        raw=read(joinpath(@__DIR__,"fixtures/cpp_D5_export.json"),String)
        key=cache_key(P,B)
        cache=joinpath(directory,"cache");mkpath(cache)
        COAV.atomic_json(joinpath(cache,key*".json"),(;
            schema="av_checkpoint_v1",completed=true,problem=COAV.problem_descriptor(P),
            backend_commit=B.commit,executable_sha256=B.executable_sha256,
            options=PerfectFormsAdapter.backend_options(P),cache_key=key,
            raw_export_json=raw,raw_export_sha256=bytes2hex(sha256(raw))))
        err=try
            enumerate_perfect_forms(P;backend=:polyhedral_common,config=B,cache_dir=cache,work_dir=directory)
            nothing
        catch caught
            caught
        end
        @test err isa BackendError
        @test err.kind==:child_exit
        E=enumerate_perfect_forms(P;backend=:polyhedral_common,config=B,cache_dir=cache,work_dir=directory,force=false)
        @test E.provenance["from_cache"] && E.completion_status==:complete
        @test isempty(E.pending_facets)
    end
end
