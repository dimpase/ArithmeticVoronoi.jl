using ArithmeticVoronoi, Oscar, Test
const QUAV=ArithmeticVoronoi

@testset "exact fundamental units beyond the original five fields" begin
    historical=Dict(5=>(0,1),8=>(1,1),12=>(2,1),13=>(1,1),40=>(3,1))
    for D in (5,8,12,13,17,21,24,28,29,33,37,40,41,44,56,60)
        P=arithmetic_problem(D)
        units=QUAV.verified_units(P)
        e=units.epsilon
        @test abs(units.epsilon_exponent)==1
        @test e in P.order && abs(norm(e))==1
        @test tr(e)>0 && coeff(e,1)>0
        if haskey(historical,D)
            a,b=historical[D]
            @test e==a+b*P.generator
        end
        # Both signs, both directions and nontrivial powers must retain the
        # determinant coordinates used in the index-four SL extraction.
        for sign in 0:1,exponent in (-7,-1,0,1,6)
            c=QUAV.unit_coordinates(P,units,(-1)^sign*e^exponent)
            @test c.sign==sign && c.exponent==exponent
            @test c.mod_squares==sign+2mod(exponent,2)
        end
        metadata=field_metadata(P)
        @test metadata.order_is_maximal
        @test metadata.ordinary_class_number==(D in (40,60) ? "2" : "1")
    end
end
