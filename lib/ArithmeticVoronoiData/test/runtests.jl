using ArithmeticVoronoiData, Test

@testset "explicit recorded-results library" begin
    @test !isdefined(ArithmeticVoronoiData,:ArithmeticVoronoi)
    @test !isdefined(ArithmeticVoronoiData,:Oscar)
    for (collection,D,count) in ((:native_real,5,48),(:native_imaginary,-7,15),(:page,-7,15))
        catalogue=available_presentations(collection)
        @test length(catalogue)==count
        @test all(e->e.collection==collection && e.rank==2,catalogue)
        R=@test_logs (:info,r"Loading a recorded presentation") load_presentation(collection,D)
        @test R.origin==:recorded_library && R.collection==collection
        @test R.discriminant==D && R.group==:GL && !R.projective
        @test length(R.provenance.source_sha256)==64
        @test isfile(R.provenance.source_file)
        @test all(w->all(i->1<=abs(i)<=length(R.generator_matrices),w),R.relators)
        @test all(m->length(m)==R.rank && all(r->length(r)==R.rank,m),R.generator_matrices)
    end
    @test_throws ArgumentError available_presentations(:unknown)
    @test_throws ArgumentError load_presentation(:native_real,61)
    @test_throws ArgumentError load_presentation(:native_imaginary,-7;rank=3)
    @test_throws ArgumentError load_presentation(:native_real,5;group=:PGL)
end
