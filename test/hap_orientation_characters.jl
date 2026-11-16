using ArithmeticVoronoi, Oscar, Test

@testset "Trivial and reversing orientation characters" begin
    GAP=Oscar.GAP
    @test GAP.Packages.load(get(ENV,"ARITHMETIC_VORONOI_HAP_ROOT","hap"))
    GAP.Globals.Read(GAP.GapObj(joinpath(ArithmeticVoronoi.PACKAGE_ROOT,"gap","ArithmeticVoronoi.g")))
    I=[1 0;0 1];r=[-1 0;0 1];u=[1 1;0 1]
    wire(A)=[collect(row) for row in eachrow(A)]
    elements=wire.([I,r,u,[1 -1;0 1],r*u,u*r])
    # Two orientation modules for the same reflection stabilizer. Transporters
    # include elements outside it, where the coset correction matters.
    data=Dict("dimensions"=>[1,0,2],"elts"=>elements,
        "stabilizers"=>[[[wire(r)]],[],[[wire(r)],[wire(r)]]],
        "evenStabilizers"=>[[[wire(r)]],[],[[wire(r)],[wire(I)]]],
        "boundaries"=>[[[]],[],[[],[]]],"generators"=>wire.([r,u]),
        "name"=>"orientation character test","discriminant"=>1,"rank"=>2,"groupKind"=>"GL")
    K=GAP.Globals.HAPArithmeticVoronoiComplex(GAP.GapObj(data;recursive=true),2)
    action=GAP.evalstr("K -> K!.action")(K)
    original=GAP.evalstr("""
    function(K,E,k,g)
        local H,r,u;
        H:=K!.stabilizer(2,k);
        r:=CanonicalRightCosetElement(H,One(H))^-1
            *CanonicalRightCosetElement(H,K!.elts[g]^-1);
        u:=r*K!.elts[g];
        if u in E then return 1; fi;
        return -1;
    end
    """)
    even=[GAP.Globals.Group(GAP.GapObj([wire(r)];recursive=true)),
          GAP.Globals.Group(GAP.GapObj([wire(I)];recursive=true))]
    for g in eachindex(elements)
        @test Int(action(2,1,g))==1
        for k in 1:2;@test action(2,k,g)==original(K,even[k],k,g);end
    end
    @test Int(action(2,2,2))==-1
end
