using ArithmeticVoronoi, JSON3, Test
include(joinpath(@__DIR__,"../scripts/request_arguments.jl"))

@testset "arithmetic request records" begin
    q=voronoi_request(-7)
    @test q["field"]["discriminant"]==-7
    @test normalize_voronoi_request(JSON3.read(JSON3.write(q)))==q
    @test normalize_voronoi_request((;group=:GL,rank=2,
        field=(;kind="quadratic",radicand=-7),order=(;kind="maximal"),lattice=(;kind="free")))==q
    for (d,D) in ((-1,-4),(2,8),(3,12),(5,5))
        r=deepcopy(q);r["field"]=Dict("kind"=>"quadratic","radicand"=>d)
        @test normalize_voronoi_request(r)==voronoi_request(D)
    end
    for D in (0,1,-12,20,45)
        @test_throws ArgumentError voronoi_request(D)
    end
    for (key,value) in (("rank",1),("rank",2.0),("rank",true),("group","PSL"),
        ("field",Dict("kind"=>"number_field","polynomial"=>[-1,-1,0,1])),
        ("order",Dict("kind"=>"nonmaximal")),("lattice",Dict("kind"=>"projective")),
        ("lattice",Dict("kind"=>"free","ideal"=>"nonprincipal")),
        ("field",Dict("kind"=>"quadratic","radicand"=>-12)),
        ("field",Dict("kind"=>"quadratic","radicand"=>-7,"discriminant"=>-7)),
        ("field",Dict("kind"=>"quadratic")))
        r=deepcopy(q);r[key]=value
        @test_throws ArgumentError normalize_voronoi_request(r)
    end
    r=deepcopy(q);r["subgroup"]="Gamma0"
    @test_throws ArgumentError normalize_voronoi_request(r)
    r=deepcopy(q);delete!(r,"order")
    @test_throws ArgumentError normalize_voronoi_request(r)
    @test_throws ArgumentError enumerate_perfect_forms(voronoi_request(-7;group=:SL))
    @test_throws ArgumentError hap_voronoi_complex(voronoi_request(5;group=:SL))
    @test_throws ArgumentError hap_voronoi_complex(q;resolution_length=0)
    withenv("PERFECT_FORMS_BACKEND_CONFIG"=>"/no/such/config.json") do
        E=enumerate_perfect_forms(q)
        @test E.problem.discriminant==-7
        @test E.provenance["backend"]=="native"
        @test length(E.representatives)==1
    end
    mktempdir() do directory
        path=joinpath(directory,"requests.json")
        write(path,JSON3.write([q,voronoi_request(5)]))
        options=request_arguments(["--requests=$path","-7","--resume"];default_discriminants=[8],default_output=directory)
        @test options.requests==[q,voronoi_request(5)]
        @test options.resume
        @test matching_request(Dict("request"=>q),q)
        @test !matching_request(Dict("request"=>voronoi_request(5)),q)
        @test !matching_request(Dict(),q)
        @test request_arguments(String[];default_discriminants=[8],default_output=directory).requests==[voronoi_request(8)]
        write(path,"[]")
        @test_throws ArgumentError request_arguments(["--requests=$path"];default_discriminants=[8],default_output=directory)
        write(path,JSON3.write(q))
        @test_throws ArgumentError request_arguments(["--requests=$path"];default_discriminants=[8],default_output=directory)
        r=deepcopy(q);r["rank"]=1
        output=joinpath(directory,"invalid")
        @test_throws ArgumentError run_field_pipeline(r;directory=output)
        @test !ispath(output)
    end
end
