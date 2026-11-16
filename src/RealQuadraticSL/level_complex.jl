# Every simplex is a flag of Voronoi cells, ordered by dimension.
# Its level section is a product of simplices. Flag stabilizers fix it pointwise.
struct FlagOrbit
    flag::Vector{MorseVertex}
    stabilizer::Vector{Any}
end
struct LevelCell
    flag_orbit::Int
    translate::Any
    exponent::Int
    sign::Int
end
struct MorseLevelComplex
    experiment::MorseExperiment
    certificate::Any
    flag_orbits::Vector{FlagOrbit}
    cells::Vector{Vector{LevelCell}}
    stabilizers::Vector{Vector{Vector{Any}}}
    boundaries::Vector{Vector{Vector{Tuple{Int,Any}}}}
end

function all_flag_orbits(M)
    C=M.complex
    orbits=FlagOrbit[]
    lookup=[[Dict{Any,Tuple{Int,Any}}() for c in layer] for layer in C.cells]
    for (k,layer) in enumerate(M.flags),(i,fullflags) in enumerate(layer)
        table=lookup[k][i]
        for full in fullflags,mask in 0:(1<<(length(full)-1))-1
            flag=[full[j] for j in 1:length(full)-1 if !iszero(mask & (1<<(j-1)))]
            push!(flag,last(full))
            key=Tuple(v.key for v in flag)
            haskey(table,key) && continue
            stabilizer=Any[]
            id=length(orbits)+1
            for s in C.stabilizers[k][i]
                image=Tuple(morse_act(M,s,v).key for v in flag)
                table[image]=(id,s)
                image==key && push!(stabilizer,s)
            end
            push!(orbits,FlagOrbit(flag,stabilizer))
        end
    end
    orbits,lookup
end

function build_level_complex(M,data,certificate;progress=false)
    @assert certificate.passed
    C=M.complex;P=C.enumeration.problem
    ids=[[sum(length,C.cells[1:k-1];init=0)+i for i in eachindex(layer)] for (k,layer) in enumerate(C.cells)]
    height(v)=certificate.bases[ids[v.dim+1][v.orbit]]+v.shift
    orbits,lookup=all_flag_orbits(M)
    if progress;println("Barycentric simplex orbits: ",length(orbits));flush(stdout);end
    top=length(C.cells)-2
    cells=[LevelCell[] for _ in 0:top]
    stabilizers=[Vector{Any}[] for _ in 0:top]
    cell_ids=Dict{Tuple{Int,Int,Int},Int}()
    flip=[any(s->det(s)==-1,o.stabilizer) for o in orbits]
    for (i,o) in enumerate(orbits)
        length(o.flag)>=2 || continue
        h=height.(o.flag)
        @assert all(!isinteger,h)
        lo=floor(Int,-maximum(h))+1
        hi=ceil(Int,-minimum(h))-1
        k=length(o.flag)-1
        for q in lo:hi,sgn in (flip[i] ? (1,) : (1,-1))
            t=diagonal_matrix(P.field,[sgn*M.units.epsilon^q;fill(one(P.field),P.rank-1)])
            push!(cells[k],LevelCell(i,t,q,sgn))
            push!(stabilizers[k],[t*s*inv(t) for s in o.stabilizer if det(s)==1])
            cell_ids[(i,q,sgn)]=length(cells[k])
        end
    end
    if progress;println("SL level cell orbits: ",length.(cells));flush(stdout);end
    boundaries=[[Tuple{Int,Any}[] for c in layer] for layer in cells]
    for k in 2:length(cells),(i,c) in enumerate(cells[k])
        F=orbits[c.flag_orbit].flag
        for j in eachindex(F)
            sub=[F[l] for l in eachindex(F) if l!=j]
            h=height.(sub).+c.exponent
            minimum(h)<0<maximum(h) || continue
            v=last(sub)
            key=Tuple(morse_act(M,inv(v.translate),w).key for w in sub)
            orbit,s=lookup[v.dim+1][v.orbit][key]
            t=c.translate*v.translate*s
            u=AV.unit_coordinates(P,M.units,det(t))
            sgn=flip[orbit] ? 1 : (-1)^u.sign
            number=cell_ids[(orbit,Int(u.exponent),sgn)]
            target=cells[k-1][number]
            correction=first(x for x in orbits[orbit].stabilizer if det(t)*det(x)==det(target.translate))
            g=t*correction*inv(target.translate)
            @assert det(g)==1
            # ∂(level(simplex)) = -level(∂simplex), with the height coorientation.
            # Zero-dimensional cells use the positive augmentation orientation.
            lower_orientation=length(sub)==2 ? Int(sign(height(sub[2])-height(sub[1]))) : 1
            push!(boundaries[k][i],((-1)^j*lower_orientation*number,g))
        end
    end
    MorseLevelComplex(M,certificate,orbits,cells,stabilizers,boundaries)
end

function validate_level_complex(L)
    P=L.experiment.complex.enumeration.problem
    I=identity_matrix(P.field,P.rank)
    @assert all(layer->all(S->I in S && all(g->det(g)==1,S),layer),L.stabilizers)
    function term_key(k,i,g)
        minimum(AV.matrix_key(AV.encode_matrix(P,g*s)) for s in L.stabilizers[k][i])
    end
    for k in 2:length(L.cells),(i,c) in enumerate(L.cells[k])
        b=L.boundaries[k][i]
        @assert !isempty(b)
        @assert all(t->det(t[2])==1,b)
        if k==2
            @assert sum(sign(j) for (j,g) in b)==0
        else
            coefficients=Dict{Any,Int}()
            for (j,g) in b,(l,h) in L.boundaries[k-1][abs(j)]
                key=(abs(l),term_key(k-2,abs(l),g*h))
                coefficients[key]=get(coefficients,key,0)+sign(j)*sign(l)
            end
            @assert all(iszero,values(coefficients)) "equivariant boundary squared is nonzero"
        end
        expected=Dict{Any,Int}()
        for (j,g) in b
            key=(abs(j),term_key(k-1,abs(j),g))
            expected[key]=get(expected,key,0)+sign(j)
        end
        for s in L.stabilizers[k][i]
            actual=Dict{Any,Int}()
            for (j,g) in b
                key=(abs(j),term_key(k-1,abs(j),s*g))
                actual[key]=get(actual,key,0)+sign(j)
            end
            @assert actual==expected "boundary is not invariant under its stabilizer"
        end
    end
    true
end

function level_hap_complex(L;resolution_length=4,hap_root=get(ENV,"ARITHMETIC_VORONOI_HAP_ROOT","hap"))
    @assert validate_level_complex(L)
    P=L.experiment.complex.enumeration.problem
    GAP=Oscar.GAP
    GAP.Packages.load(hap_root) || error("HAP not available")
    GAP.Globals.Read(GAP.GapObj(joinpath(AV.PACKAGE_ROOT,"gap","ArithmeticVoronoi.g")))
    elts=Any[identity_matrix(P.field,P.rank)]
    indices=Dict(AV.matrix_key(AV.encode_matrix(P,only(elts)))=>1)
    function index(g)
        key=AV.matrix_key(AV.encode_matrix(P,g))
        get!(indices,key) do
            push!(elts,g);length(elts)
        end
    end
    for layer in L.stabilizers,S in layer,g in S;index(g);end
    boundaries=[[[[j,index(g)] for (j,g) in b] for b in layer] for layer in L.boundaries]
    stabs=[[AV.cell_generators(P,S) for S in layer] for layer in L.stabilizers]
    generators=hasproperty(L,:generators) ? L.generators : unique([g for layer in stabs for S in layer for g in S] ∪ [g for layer in L.boundaries for b in layer for (j,g) in b])
    wire(g)=[[BigInt(ZZ(x)) for x in row] for row in AV.matrix_rows(AV.encode_matrix(P,g))]
    data=Dict("dimensions"=>length.(L.cells),"elts"=>wire.(elts),
        "stabilizers"=>[[wire.(S) for S in layer] for layer in stabs],
        "evenStabilizers"=>[[wire.(S) for S in layer] for layer in stabs],
        "boundaries"=>boundaries,"generators"=>wire.(generators),
        "name"=>"SL($(P.rank),O), D=$(P.discriminant), certified determinant-height level complex",
        "discriminant"=>P.discriminant,"rank"=>P.rank,"groupKind"=>"SL")
    GAP.Globals.HAPArithmeticVoronoiComplex(GAP.GapObj(data;recursive=true),resolution_length)
end
