# Optional comparison with VMH's original geometric orientation calculation.
# Input is a prepared Magma reference directory; nothing is loaded by the engine.
include("prepare_real_sl_magma.jl")

function compare_real_sl_magma_tangents(reference,output)
    reference,output=abspath(reference),abspath(output)
    source=joinpath(reference,"Magma")
    functions=read(joinpath(source,"Homology","Functions.m"),String)
    pattern=r"intrinsic OrientationSignByDeterminant\(.*?end intrinsic;"s
    original=match(pattern,functions)
    isnothing(original) && error("VMH orientation function not found")
    occursin("BoundaryEmbeddings",original.match) || error("reference must retain VMH's geometric orientation calculation")
    ispath(output) && error("refusing to overwrite existing test output")
    mkpath(output)
    patched=joinpath(output,"patched_functions.m")
    write(patched,functions)
    linear_vmh_tangents!(patched)
    alternative=replace(match(pattern,read(patched,String)).match,
        "OrientationSignByDeterminant"=>"KernelOrientationSign")
    literal(path)=repr(path)
    write(joinpath(output,"kernel_orientation.m"),
        "import $(literal(joinpath(source,"Initialize.m"))): dimsym;\n"*alternative*"\n")
    write(joinpath(output,"compare.m"),"""
    AttachSpec($(literal(joinpath(source,"VMH-spec"))));
    Attach("kernel_orientation.m");
    procedure RunComparison()
    SetSeed(1);
    V:=VoronoiAlgorithm();
    MinimalClasses(V);
    checks:=0;
    for layer in V`MinimalClasses do
     for F in layer do
      for g in StabilizerOfMinimalClass(F) do
       assert OrientationSignByDeterminant(F,g,V) eq KernelOrientationSign(F,g,V);
       checks+:=1;
      end for;
     end for;
    end for;
    assert checks gt 0;
    printf "Orientation comparisons passed: %o\\n",checks;
    Write("complete.txt",IntegerToString(checks));
    end procedure;
    RunComparison();
    quit;
    """)
    write(joinpath(output,"source.json"),JSON.json(Dict(
        "reference"=>reference,
        "original_functions_sha256"=>bytes2hex(sha256(functions)),
        "patched_functions_sha256"=>bytes2hex(sha256(read(patched))));pretty=4))
    open(joinpath(output,"magma.log"),"w") do io
        run(pipeline(Cmd(`magma -b compare.m`;dir=output);stdout=io,stderr=io))
    end
    complete=joinpath(output,"complete.txt")
    isfile(complete) || error("Magma comparison failed; inspect $output/magma.log")
    parse(Int,strip(read(complete,String)))
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("provide an original prepared Magma reference and a new output directory")
    println("Orientation comparisons passed: ",compare_real_sl_magma_tangents(ARGS...))
end
