# Run only through engine_isolation.jl, in a fresh Julia process.
using Oscar, JSON3, Test
const ISOLATED_ROOT=only(ARGS)
include(joinpath(ISOLATED_ROOT,"src/ArithmeticVoronoi.jl"))
using .ArithmeticVoronoi
const IAV=ArithmeticVoronoi

@testset "fresh computation and explicit checkpoint replay" begin
    @test realpath(IAV.PACKAGE_ROOT)==realpath(ISOLATED_ROOT)
    @test all(!ispath(joinpath(ISOLATED_ROOT,p)) for p in ("results","test","references","validation","lib"))
    @test !isdefined(IAV,:compare_kw)
    @test !isdefined(IAV,:kw_reference)
    @test !isdefined(IAV,:ArithmeticVoronoiData)
    withenv("PERFECT_FORMS_BACKEND_CONFIG"=>"/no/such/config.json") do
        directory=joinpath(ISOLATED_ROOT,"fresh_GL3")
        P=classical_problem(;rank=3)
        R=run_gl_pipeline(P;directory,progress=false)
        @test verify_presentation(R.gl)
        @test R.run["presentation"]["abelianization"].torsion==["2"]
        @test !R.run["resume_requested"] && R.run["replayed_loop_certificates"]==0
        ledger=JSON3.read(read(joinpath(directory,"two_cells/coverage.json"),String))
        loop=joinpath(directory,"two_cells",first(first(ledger.entries).certificates))
        bad=JSON3.read(read(loop,String),Dict{String,Any});bad["source_hash"]="corrupt"
        IAV.atomic_json(loop,bad)
        @test_throws ArithmeticValidationError run_gl_pipeline(P;directory,progress=false,resume=true)
        # A fresh computation must ignore stale/corrupt saved answers.
        fresh=run_gl_pipeline(P;directory,progress=false)
        @test fresh.run["replayed_loop_certificates"]==0
        replay=run_gl_pipeline(P;directory,progress=false,resume=true)
        @test replay.run["resume_requested"]
        @test replay.run["replayed_loop_certificates"]==replay.run["geometry"]["oriented_loops"]>0
        @test verify_presentation(replay.gl)

        # Former automatic baseline location: deliberately unusable data.
        baseline=joinpath(ISOLATED_ROOT,"results/two_cells/D5")
        mkpath(baseline);write(joinpath(baseline,"coverage.json"),"this is not JSON")
        fielddir=joinpath(ISOLATED_ROOT,"fresh_D5");mkpath(joinpath(fielddir,"two_cells"))
        write(joinpath(fielddir,"two_cells/coverage.json"),"this is not JSON")
        write(joinpath(fielddir,"GL_simplified.json"),"this is not JSON")
        F=run_field_pipeline(voronoi_request(5);directory=fielddir)
        @test F.run["backend"]=="native" && !F.enumeration.provenance["from_cache"]
        @test !F.run["resume_requested"] && F.run["replayed_loop_certificates"]==0
        @test F.run["presentations"]["GL"]["abelianization"].torsion==["2"]
        @test F.run["presentations"]["GL"]["abelianization"].free_rank==1
        @test read(joinpath(baseline,"coverage.json"),String)=="this is not JSON"
        @test verify_presentation(F.gl) && verify_sl_presentation(F.sl)
        g=matrix(F.enumeration.problem.field,[1 1;0 1])
        @test verify_word_certificate(F,g,matrix_to_word(F,g;target=:GL))
        badsimple=JSON3.read(read(joinpath(fielddir,"GL_simplified.json"),String),Dict{String,Any})
        badsimple["old_generator_count"]+=1
        IAV.atomic_json(joinpath(fielddir,"GL_simplified.json"),badsimple)
        resumed=run_field_pipeline(voronoi_request(5);directory=fielddir,resume=true)
        @test resumed.run["replayed_loop_certificates"]==2sum(length,resumed.context.ridges)
        saved=JSON3.read(read(joinpath(fielddir,"GL_simplified.json"),String))
        @test saved.old_generator_count==length(resumed.gl.generator_matrices)
        @test verify_presentation(resumed.gl) && verify_sl_presentation(resumed.sl)
    end
end
