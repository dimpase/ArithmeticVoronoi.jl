# Optional independent references; never called by the engine or package loading.
using JSON
length(ARGS)==2 || error("provide a VMH-TotallyRealNumberFields checkout and a durable output directory")
source,output=abspath.(ARGS)
isfile(joinpath(source,"Magma","VMH-spec")) || error("VMH source directory not found")
revision=readchomp(`git -C $source rev-parse HEAD`)
for (D,d) in ((5,5),(12,3))
    target=joinpath(output,"magma_D$D")
    ispath(target) && error("refusing to append Magma output to existing $target")
    mkpath(target)
    for directory in ("Magma","GAP")
        cp(joinpath(source,directory),joinpath(target,directory))
    end
    write(joinpath(target,"Magma","BasicData.m"),
        "K:=QuadraticField($d); n:=2; V:=KMatrixSpace(K,n,1); L:=[IntegralBasis(K)[i]*V.j: i in [1..Degree(K)],j in [1..n]];\n")
    write(joinpath(target,"compute.m"),"""
    AttachSpec("Magma/VMH-spec");
    SetSeed(1);
    V:=VoronoiAlgorithm();
    ComputeComplexGL("complex_gl.g",V);
    Write("complete.txt","complete");
    quit;
    """)
    write(joinpath(target,"source.json"),JSON.json(Dict("repository"=>source,"revision"=>revision,"discriminant"=>D)))
    open(joinpath(target,"magma.log"),"w") do io
        run(pipeline(Cmd(`magma -b compute.m`;dir=target);stdout=io,stderr=io))
    end
    isfile(joinpath(target,"complete.txt")) || error("Magma did not finish; inspect $target/magma.log")
end
