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

# A reference for H_k needs cells only through degree k+1. Truncating the
# two dimension loops leaves all requested cells and boundary computations
# unchanged; this option does not assert contractibility of the skeleton.
function truncate_vmh_dimensions!(path,dimension)
    text=read(path,String)
    for condition in ("dimsym-codim ge 1", "dimsym-codim ge 0")
        original="while $condition do"
        count(original,text)==1 || error("VMH MinClass.m changed; review dimension truncation")
        text=replace(text,original=>"while $condition and codim le $dimension do")
    end
    write(path,text)
end

# Cache only values independent of the acting matrix. These are the same
# VMH stabilizer and affine tangent-space computations, evaluated once per
# minimal class; no native Julia geometric data enter the reference.
function memoize_vmh_cells!(directory)
    types=joinpath(directory,"Types.m")
    text=read(types,String)
    old="Aut,CanonicalFormMinimalClass;"
    count(old,text)==1 || error("VMH Types.m changed; review memoization")
    write(types,replace(text,old=>"Aut,CanonicalFormMinimalClass,MinimalClassAut,OrientationTangent;"))
    classes=joinpath(directory,"ReQuad","Voronoi","MinClass.m")
    text=read(classes,String)
    old=" return AutomorphismGroup(NewDAForm((CanonicalFormOfMinimalClass(A)`Matrix)^-1));"
    count(old,text)==1 || error("VMH stabilizer source changed; review memoization")
    write(classes,replace(text,old=>" if not assigned(A`MinimalClassAut) then\n  A`MinimalClassAut:=AutomorphismGroup(NewDAForm((CanonicalFormOfMinimalClass(A)`Matrix)^-1));\n end if;\n return A`MinimalClassAut;"))
    functions=joinpath(directory,"Homology","Functions.m")
    text=read(functions,String)
    first_line=" if not assigned(W`MinimalClasses) then"
    last_line=" Images:=[Dagger(g)*GensOfTSpace[i]*(g): i in Indices];"
    count(first_line,text)==1 && count(last_line,text)==1 || error("VMH orientation source changed; review memoization")
    text=replace(text,first_line=>" if not assigned(F`OrientationTangent) then\n"*first_line,
        last_line=>" F`OrientationTangent:=<Space,[GensOfTSpace[i]:i in Indices]>;\n end if;\n Space:=F`OrientationTangent[1];\n Images:=[Dagger(g)*basis*g: basis in F`OrientationTangent[2]];")
    write(functions,text)
end

function prepare_real_sl_magma(source,output;cases=((D=5,rank=2),(D=12,rank=2)),
        deduplicate_faces=false,max_cell_dimension=nothing,memoize_cell_data=false)
    isnothing(max_cell_dimension) ||
        (max_cell_dimension isa Integer && !(max_cell_dimension isa Bool) && max_cell_dimension>=2) ||
        throw(ArgumentError("max_cell_dimension must be an integer at least two"))
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
        source_files=("Types.m","Homology/Functions.m","ReQuad/Voronoi/MinClass.m")
        source_hashes=Dict(file=>bytes2hex(sha256(read(joinpath(target,"Magma",file)))) for file in source_files)
        original_hash=bytes2hex(sha256(read(class_source)))
        deduplicate_faces && deduplicate_vmh_faces!(class_source)
        memoize_cell_data && memoize_vmh_cells!(joinpath(target,"Magma"))
        isnothing(max_cell_dimension) || truncate_vmh_dimensions!(class_source,max_cell_dimension)
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
            "max_cell_dimension"=>max_cell_dimension,"exact_face_deduplication"=>deduplicate_faces,
            "memoized_cell_data"=>memoize_cell_data,"original_source_sha256"=>source_hashes,
            "used_source_sha256"=>Dict(file=>bytes2hex(sha256(read(joinpath(target,"Magma",file)))) for file in source_files),"original_MinClass_sha256"=>original_hash,
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
