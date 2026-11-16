# Optional independent references; never called by the engine or package loading.
using ArithmeticVoronoi, JSON, SHA

# Optional performance-only change in the copied reference source. A face is
# determined by its intersection of minimum-vector sets: repeated facet tuples
# with the same intersection need not be expanded again. Equivalence testing,
# stabilizers, orientations and the exported complex remain VMH computations.
function deduplicate_vmh_faces!(path)
    text=read(path,String)
    replacements=("FaceTuples:=[**];"=>"FaceTuples:=[**]; SeenIntersections:={};",
        "if #Intersection gt 1 then"=>"if #Intersection gt 1 and Intersection notin SeenIntersections then\n     Include(~SeenIntersections,Intersection);",
        "if #IntersectionTemp gt 1 then"=>"if #IntersectionTemp gt 1 and IntersectionTemp notin SeenIntersections then\n      Include(~SeenIntersections,IntersectionTemp);")
    for (replacement,expected) in zip(replacements,(2,1,1))
        count(first(replacement),text)==expected || error("VMH MinClass.m changed; review the optional deduplication patch")
        text=replace(text,replacement)
    end
    write(path,text)
end

function prepare_real_sl_magma(source,output;cases=((D=5,rank=2),(D=12,rank=2)),deduplicate_faces=false)
    source,output=abspath(source),abspath(output)
    isfile(joinpath(source,"Magma","VMH-spec")) || error("VMH source directory not found")
    revision=readchomp(`git -C $source rev-parse HEAD`)
    for case in cases
        P=arithmetic_problem(case.D;rank=case.rank)
        D=P.discriminant;n=P.rank;d=mod(D,4)==1 ? D : div(D,4)
        target=joinpath(output,n==2 ? "magma_D$D" : "magma_rank$(n)_D$D")
        ispath(target) && error("refusing to append Magma output to existing $target")
        mkpath(target)
        for directory in ("Magma","GAP")
            cp(joinpath(source,directory),joinpath(target,directory))
        end
        class_source=joinpath(target,"Magma","ReQuad","Voronoi","MinClass.m")
        original_hash=bytes2hex(sha256(read(class_source)))
        deduplicate_faces && deduplicate_vmh_faces!(class_source)
        write(joinpath(target,"Magma","BasicData.m"),
            "K:=QuadraticField($d); n:=$n; V:=KMatrixSpace(K,n,1); L:=[IntegralBasis(K)[i]*V.j: i in [1..Degree(K)],j in [1..n]];\n")
        write(joinpath(target,"compute.m"),"""
        AttachSpec("Magma/VMH-spec");
        SetSeed(1);
        V:=VoronoiAlgorithm();
        ComputeComplexGL("complex_gl.g",V);
        Write("complete.txt","complete");
        quit;
        """)
        write(joinpath(target,"source.json"),JSON.json(Dict("repository"=>source,"revision"=>revision,"discriminant"=>D,"rank"=>n,
            "exact_face_deduplication"=>deduplicate_faces,"original_MinClass_sha256"=>original_hash,
            "used_MinClass_sha256"=>bytes2hex(sha256(read(class_source))))))
        open(joinpath(target,"magma.log"),"w") do io
            run(pipeline(Cmd(`magma -b compute.m`;dir=target);stdout=io,stderr=io))
        end
        isfile(joinpath(target,"complete.txt")) || error("Magma did not finish; inspect $target/magma.log")
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("provide a VMH-TotallyRealNumberFields checkout and a durable output directory")
    prepare_real_sl_magma(ARGS...)
end
