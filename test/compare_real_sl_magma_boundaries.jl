# Optional exact comparison with VMH's original exhaustive boundary search.
# Every vector and minimal class is computed by Magma, independently of Julia.
include("prepare_real_sl_magma.jl")

function compare_real_sl_magma_boundaries(reference,output;configurations_only=false)
    reference,output=abspath(reference),abspath(output)
    source=joinpath(reference,"Magma")
    functions=read(joinpath(source,"Homology","Functions.m"),String)
    pattern=r"intrinsic BoundaryEmbeddings\(.*?end intrinsic;"s
    original=match(pattern,functions)
    isnothing(original) && error("VMH boundary function not found")
    occursin("ImageSets:=Subsets",original.match) || error("reference must retain VMH's original ordered-basis search")
    ispath(output) && error("refusing to overwrite existing comparison output")
    mkpath(output)
    patched=joinpath(output,"patched_functions.m")
    write(patched,functions)
    incremental_vmh_boundaries!(patched)
    alternative=replace(match(pattern,read(patched,String)).match,
        "BoundaryEmbeddings"=>"IncrementalBoundaryEmbeddings")
    literal(path)=repr(path)
    implementation="import $(literal(joinpath(source,"Initialize.m"))): n;\n"*alternative*"\n"
    write(joinpath(output,"incremental_boundary.m"),implementation)
    setup=configurations_only ? """
    // Small positive forms exercise arbitrary-rank linear dependencies without
    // using a stored Voronoi enumeration or a known boundary answer.
    I:=IdentityMatrix(K,n);
    Block:=I; Block[1,2]:=-1/2; Block[2,1]:=-1/2;
    Path:=I;
    for i in [1..n-1] do Path[i,i+1]:=-1/2; Path[i+1,i]:=-1/2; end for;
    Dense:=(I+Matrix(K,n,n,[1:i in [1..n*n]]))/2;
    Forms:=[NewDAForm(A):A in [I,Block,Path,Dense]];
    Pairs:=[<F,G>:F in Forms,G in Forms];
    """ : """
    V:=VoronoiAlgorithm();
    MinimalClasses(V);
    Pairs:=[];
    for d in [2..#V`MinimalClasses] do
     for F in V`MinimalClasses[d] do
      for G in V`MinimalClasses[d-1] do Append(~Pairs,<F,G>); end for;
     end for;
    end for;
    """
    write(joinpath(output,"compare.m"),"""
    import $(literal(joinpath(source,"BasicData.m"))): K;
    import $(literal(joinpath(source,"Initialize.m"))): n;
    AttachSpec($(literal(joinpath(source,"VMH-spec"))));
    Attach("incremental_boundary.m");
    procedure RunComparison()
    SetSeed(1);
    $setup
    checks:=0; oldtime:=0.0; newtime:=0.0;
    for pair in Pairs do
     Large:=pair[1]; Small:=pair[2];
     t:=Cputime(); old:=BoundaryEmbeddings(Large,Small,IsIntegralUnit); oldtime+:=Cputime(t);
     t:=Cputime(); new:=IncrementalBoundaryEmbeddings(Large,Small,IsIntegralUnit); newtime+:=Cputime(t);
     S:=StabilizerOfMinimalClass(Small);
     assert #old eq #new;
     assert forall{g:g in old|exists{h:h in new|h^-1*g in S}};
     checks+:=1;
    end for;
    assert checks eq #Pairs and checks gt 0;
    printf "Boundary comparisons passed: %o; original %o s; incremental %o s\\n",checks,oldtime,newtime;
    Write("complete.txt",IntegerToString(checks));
    Write("timings.txt",Sprint(<oldtime,newtime>));
    end procedure;
    RunComparison();
    quit;
    """)
    write(joinpath(output,"source.json"),JSON.json(Dict(
        "reference"=>reference,"configurations_only"=>configurations_only,
        "original_functions_sha256"=>bytes2hex(sha256(functions)),
        "incremental_function_sha256"=>bytes2hex(sha256(implementation)));pretty=4))
    open(joinpath(output,"magma.log"),"w") do io
        run(pipeline(Cmd(`magma -b compare.m`;dir=output);stdout=io,stderr=io))
    end
    complete=joinpath(output,"complete.txt")
    isfile(complete) || error("Magma comparison failed; inspect $output/magma.log")
    parse(Int,strip(read(complete,String)))
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (2,3) || error("provide an original prepared reference, a new output directory, and optionally --configurations-only")
    configurations_only=length(ARGS)==3
    configurations_only && ARGS[3]!="--configurations-only" && error("unknown comparison option")
    println("Boundary comparisons passed: ",compare_real_sl_magma_boundaries(ARGS[1],ARGS[2];configurations_only))
end
