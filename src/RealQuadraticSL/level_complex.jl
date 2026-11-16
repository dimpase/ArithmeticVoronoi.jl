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

function all_flag_orbits(M;faces=morse_cell_faces(M),max_length=length(M.complex.cells),lookup_max_length=max_length,height=nothing,progress=false)
    C=M.complex;P=C.enumeration.problem
    orbits=FlagOrbit[]
    lookup=[[Dict{Any,Tuple{Int,Any}}() for c in layer] for layer in C.cells]
    I=identity_matrix(P.field,P.rank)
    crosses(lo,hi)=floor(BigInt,-hi)+1<=ceil(BigInt,-lo)-1
    meets_level(flag)=isnothing(height) || crosses(minimum(height,flag),maximum(height,flag))
    bounds=isnothing(height) ? nothing : [[begin
        lo=hi=height(M.vertices[k][i])
        for v in faces[k][i];h=height(v);lo=min(lo,h);hi=max(hi,h);end
        (lo,hi)
    end for i in eachindex(layer)] for (k,layer) in enumerate(C.cells)]
    for (k,layer) in enumerate(C.cells),i in eachindex(layer)
        table=lookup[k][i]
        function descend(flag,S)
            if meets_level(flag)
                id=length(orbits)+1
                # Boundary lookup needs only proper subflags. In a large
                # skeleton, indexing all translates of its top cells costs
                # far more space than retaining their orbit representatives.
                if length(flag)<=lookup_max_length
                    for s in C.stabilizers[k][i]
                        image=Tuple(morse_act(M,s,v).key for v in flag)
                        @assert !haskey(table,image) || first(table[image])==id "a partial-flag orbit was repeated"
                        table[image]=(id,s)
                    end
                end
                push!(orbits,FlagOrbit(flag,S))
            end
            length(flag)>=max_length && return
            v=first(flag)
            v.dim==0 && return
            if !isnothing(bounds)
                lo,hi=bounds[v.dim+1][v.orbit]
                lo=min(minimum(height,flag),lo+v.shift)
                hi=max(maximum(height,flag),hi+v.shift)
                # Every extension lies among these faces. If their entire
                # height interval misses the level, no subflag can meet it.
                crosses(lo,hi) || return
            end
            candidates=[morse_act(M,v.translate,w) for w in faces[v.dim+1][v.orbit]]
            seen=Set{Any}()
            scalar=all(s->s==I || s==-I,S)
            for w in candidates
                w.key in seen && continue
                next=Any[]
                if scalar
                    next=S
                else
                    for s in S
                        image=morse_act(M,s,w)
                        push!(seen,image.key)
                        image.key==w.key && push!(next,s)
                    end
                end
                descend([w;flag],next)
            end
        end
        descend([M.vertices[k][i]],C.stabilizers[k][i])
        if progress
            println("Partial flag orbits through cell ",(k-1,i),": ",length(orbits));flush(stdout)
        end
        # Release temporary matrices and action lookups. Interned vertices
        # remain shared by the incidence posets and flag representatives.
        clear_action_cache!(M)
    end
    orbits,lookup
end

function build_level_complex(M,data,certificate;progress=false,max_cell_dimension=nothing,faces=morse_cell_faces(M))
    @assert certificate.passed
    C=M.complex;P=C.enumeration.problem
    ids=[[sum(length,C.cells[1:k-1];init=0)+i for i in eachindex(layer)] for (k,layer) in enumerate(C.cells)]
    heights=Dict{NTuple{3,Int},Rational{BigInt}}()
    height(v)=get!(heights,(v.dim,v.orbit,v.shift)) do
        certificate.bases[ids[v.dim+1][v.orbit]]+v.shift
    end
    top=length(C.cells)-2
    if !isnothing(max_cell_dimension)
        max_cell_dimension isa Integer && !(max_cell_dimension isa Bool) && max_cell_dimension>=0 ||
            throw(ArgumentError("max_cell_dimension must be a nonnegative integer"))
        top=min(top,max_cell_dimension)
    end
    orbits,lookup=all_flag_orbits(M;faces,max_length=top+2,lookup_max_length=top+1,height,progress)
    if progress;println("Level-intersecting barycentric simplex orbits: ",length(orbits));flush(stdout);end
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

# Flag stabilizers fix every vertex. Coarser cell models can override this
# character when a stabilizer reverses the orientation of its cell.
level_orientation(L::MorseLevelComplex,k,i,g)=1

function validate_level_complex(L)
    P=L.experiment.complex.enumeration.problem;I=identity_matrix(P.field,P.rank)
    centre=iseven(P.rank) ? (I,-I) : (I,)
    for k in eachindex(L.cells),i in eachindex(L.cells[k])
        S=L.stabilizers[k][i]
        @assert all(z->z in S,centre) && all(g->det(g)==1,S)
        k==1 && @assert all(g->level_orientation(L,k,i,g)==1,S)
    end
    function add!(coefficients,k,j,g,value)
        products=[g*s for s in L.stabilizers[k][abs(j)]]
        keys=[AV.matrix_key(AV.encode_matrix(P,h)) for h in products]
        chosen=argmin(keys);key=(abs(j),keys[chosen])
        factor=level_orientation(L,k,abs(j),L.stabilizers[k][abs(j)][chosen])
        coefficients[key]=get(coefficients,key,0)+sign(j)*factor*value
    end
    for k in 2:length(L.cells),i in eachindex(L.cells[k])
        b=L.boundaries[k][i]
        L isa MorseLevelComplex && @assert !isempty(b)
        @assert all(t->det(t[2])==1,b)
        if k==2;@assert sum(sign(j) for (j,g) in b)==0;end
        if k>=3
            coefficients=Dict{Any,Int}()
            for (j,g) in b,(l,h) in L.boundaries[k-1][abs(j)]
                add!(coefficients,k-2,l,g*h,sign(j))
            end
            @assert all(iszero,values(coefficients)) "equivariant boundary squared is nonzero"
        end
        expected=Dict{Any,Int}()
        for (j,g) in b;add!(expected,k-1,j,g,1);end
        filter!(p->last(p)!=0,expected)
        for s in L.stabilizers[k][i]
            actual=Dict{Any,Int}()
            for (j,g) in b;add!(actual,k-1,j,s*g,1);end
            filter!(p->last(p)!=0,actual)
            factor=level_orientation(L,k,i,s)
            @assert actual==Dict(key=>factor*v for (key,v) in expected) "equivariant boundary has incorrect stabilizer action"
        end
    end
    true
end

function level_hap_complex(L;resolution_length=4,hap_root=get(ENV,"ARITHMETIC_VORONOI_HAP_ROOT","hap"))
    @assert validate_level_complex(L)
    computed_dimension=length(L.cells)-1
    certified_dimension=length(L.experiment.complex.cells)-2
    if computed_dimension<certified_dimension && resolution_length>computed_dimension
        throw(ArgumentError("the computed skeleton supports resolutions only through degree $computed_dimension"))
    end
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
    evens=[[AV.cell_generators(P,[g for g in S if level_orientation(L,k,i,g)==1])
            for (i,S) in enumerate(layer)] for (k,layer) in enumerate(L.stabilizers)]
    generators=hasproperty(L,:generators) ? L.generators : unique([g for layer in stabs for S in layer for g in S] ∪ [g for layer in L.boundaries for b in layer for (j,g) in b])
    wire(g)=[[BigInt(ZZ(x)) for x in row] for row in AV.matrix_rows(AV.encode_matrix(P,g))]
    data=Dict("dimensions"=>length.(L.cells),"elts"=>wire.(elts),
        "stabilizers"=>[[wire.(S) for S in layer] for layer in stabs],
        "evenStabilizers"=>[[wire.(S) for S in layer] for layer in evens],
        "boundaries"=>boundaries,"generators"=>wire.(generators),
        "name"=>"SL($(P.rank),O), D=$(P.discriminant), certified determinant-height level complex",
        "discriminant"=>P.discriminant,"rank"=>P.rank,"groupKind"=>"SL")
    data["certifiedCellDimension"]=certified_dimension
    data["computedCellDimension"]=computed_dimension
    GAP.Globals.HAPArithmeticVoronoiComplex(GAP.GapObj(data;recursive=true),resolution_length)
end
