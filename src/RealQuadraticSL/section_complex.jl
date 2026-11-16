# Certified cell sections and compatible local subdivision of the regular level.
struct SectionOrbit
    parent_dimension::Int
    parent_orbit::Int
    exponent::Int
    component::Int
    witness::MorseVertex
    stabilizer::Vector{Any}
    augmentation::Int
end
function section_components(selected,edges)
    graph=flag_adjacency(selected,edges)
    component=zeros(Int,isempty(selected) ? 0 : maximum(selected))
    members=Vector{Int}[]
    for v in sort!(collect(selected))
        component[v]>0 && continue
        push!(members,[v]);c=length(members);component[v]=c
        for u in members[c],w in graph[u]
            component[w]>0 && continue
            component[w]=c;push!(members[c],w)
        end
    end
    members,component
end

function section_symmetries(M,vertices,stabilizer)
    P=M.complex.enumeration.problem;I=identity_matrix(P.field,P.rank)
    generators=AV.cell_generators(P,stabilizer)
    filter!(g->g!=I && g!=-I,generators)
    positions=Dict(v.key=>j for (j,v) in enumerate(vertices))
    permutations=Vector{Int}[]
    for g in generators
        A=AV.encode_matrix(P,g)
        rows=Tuple(Tuple(compact_integer(BigInt(ZZ(A[a,b]))) for b in 1:ncols(A)) for a in 1:nrows(A))
        small=all(row->all(x->x isa Int,row),rows);rays=Dict{Any,Any}()
        permutation=[positions[(v.dim,Tuple(sort!([get!(rays,ray) do
            morse_ray_image(rows,small,ray)
        end for ray in v.key[2]])))] for v in vertices]
        push!(permutations,permutation)
    end
    permutations
end

struct SectionLevelCell
    kind::Symbol
    orbit::Int
    translate::Any
    exponent::Int
    sign::Int
end
struct MorseSectionLevelComplex
    experiment::MorseExperiment
    certificate::Any
    section_certificates::Any
    sections::Vector{SectionOrbit}
    flag_orbits::Vector{FlagOrbit}
    cells::Vector{Vector{SectionLevelCell}}
    stabilizers::Vector{Vector{Vector{Any}}}
    boundaries::Vector{Vector{Vector{Tuple{Int,Any}}}}
    lookup::Any
    refined::Any
    parent_dimension::Int
end

# Orient maximal barycentric simplices by their successive oriented facets.
# This agrees with geometric barycentres without constructing their coordinates.
function maximal_flag_orientation(M,F)
    d=last(F).dim
    @assert length(F)==d+1 && [v.dim for v in F]==collect(0:d)
    d==0 && return 1
    top=last(F);@assert top.key==M.vertices[d+1][top.orbit].key
    a,g=only((a,g) for (a,g) in M.complex.boundaries[d+1][top.orbit]
        if morse_act(M,g,M.vertices[d][abs(a)]).key==F[end-1].key)
    lower=[morse_act(M,inv(g),v) for v in F[1:end-1]]
    (-1)^d*sign(a)*maximal_flag_orientation(M,lower)
end

# A section is a regular neighborhood of its opposite-sign boundary core.
# A point-collapse certifies a ball. Other cores retain their exact collapse
# witnesses and residual clique dimensions; they are never assumed to be balls.
function certify_level_sections(M,certificate;faces=morse_cell_faces(M),progress=false,on_section_certificate=nothing)
    @assert certificate.passed
    C=M.complex;P=C.enumeration.problem;I=identity_matrix(P.field,P.rank)
    ids=[[sum(length,C.cells[1:k-1];init=0)+i for i in eachindex(layer)] for (k,layer) in enumerate(C.cells)]
    height(v)=certificate.bases[ids[v.dim+1][v.orbit]]+v.shift
    lower=[Dict{Any,MorseVertex}() for layer in C.cells for _ in layer]
    graphs=MorseLinkGraphs(lower,ids,M)
    sections=SectionOrbit[];proofs=Any[];marked=Set{Tuple{Int,Int,Int}}()
    lookup=[[Dict{Int,Any}() for _ in layer] for layer in C.cells]
    for k in 2:length(C.cells),i in eachindex(C.cells[k])
        vertices=faces[k][i];h=height(M.vertices[k][i]);heights=height.(vertices)
        lo=floor(Int,-max(h,maximum(heights)))+1
        hi=ceil(Int,-min(h,minimum(heights)))-1
        lo<=hi || continue
        id=ids[k][i];lower[id]=Dict(v.key=>v for v in vertices)
        link=graphs[id];empty!(lower[id])
        positions=Dict(v.key=>j for (j,v) in enumerate(vertices))
        for q in lo:hi
            direction=-Int(sign(h+q));selected=findall(x->direction*(x+q)>0,heights)
            members,component=section_components(selected,link["edges"])
            # Pad to the full face index, including unselected vertices.
            resize!(component,length(vertices))
            chosen=Set(selected)
            for v in eachindex(vertices);v in chosen || (component[v]=0);end
            maps=Dict{Int,Tuple{Int,Any}}()
            witnesses=[vertices[minimum(S)] for S in members]
            for j in eachindex(members)
                haskey(maps,j) && continue
                witness=witnesses[j];stabilizer=Any[];images=Dict{Int,Any}()
                for s in C.stabilizers[k][i]
                    target=component[positions[morse_act(M,s,witness).key]]
                    @assert target>0
                    images[target]=s
                    target==j && push!(stabilizer,s)
                end
                generators=section_symmetries(M,vertices,stabilizer)
                result=orbit_collapse_flag_graph(members[j],link["edges"],generators,link["colours"];progress)
                residual=replay_orbit_flag_reduction(members[j],link["edges"],generators,link["colours"],result.steps)
                @assert Set(keys(residual))==Set(result.residual)
                @assert result.passed==(length(residual)==1)
                core_dimension=flag_graph_dimension(residual,link["colours"])
                result.passed || push!(marked,(k,i,q))
                augmentation=1
                if k==2
                    @assert length(members[j])==1 && witness.dim==0
                    incidence=only(sign(a) for (a,g) in C.boundaries[k][i]
                        if morse_act(M,g,M.vertices[1][abs(a)]).key==witness.key)
                    augmentation=incidence*Int(sign(height(witness)-h))
                end
                push!(sections,SectionOrbit(k-1,i,q,j,witness,stabilizer,augmentation))
                number=length(sections)
                for (target,s) in images
                    @assert !haskey(maps,target)
                    maps[target]=(number,s)
                end
                record=(parent=(k-1,i),exponent=q,component=j,vertices=length(members[j]),core_dimension=core_dimension,result=result)
                push!(proofs,record)
                isnothing(on_section_certificate) || on_section_certificate(record,certificate.bases)
            end
            @assert length(maps)==length(members)
            lookup[k][i][q]=(positions=positions,component=component,witnesses=witnesses,maps=maps,
                direction=direction)
        end
        clear_action_cache!(M)
        if progress
            println("Certified section cores through original cell ",(k-1,i),": ",length(sections));flush(stdout)
        end
    end
    (;sections,proofs,lookup,refined=marked)
end

# If a d-cell section has residual core dimension r, its boundary inclusion is
# (d-r-2)-connected by PL general position and the regular-neighborhood theorem.
# Omitting only pieces with this connectivity >= m preserves exact cellular
# chains below m. Retaining all original cells through a common dimension makes
# the retained space a genuine SL-invariant subcomplex.
function section_parent_dimension(plan,full_dimension,m)
    maximum((r.parent[1] for r in plan.proofs if r.parent[1]-r.core_dimension-2<m);init=min(full_dimension,m+1))
end

function build_section_level_complex(M,certificate;faces=morse_cell_faces(M),max_cell_dimension=nothing,
        progress=false,on_section_certificate=nothing)
    full=length(M.complex.cells)-2
    if !isnothing(max_cell_dimension)
        max_cell_dimension isa Integer && !(max_cell_dimension isa Bool) && max_cell_dimension>=0 ||
            throw(ArgumentError("max_cell_dimension must be a nonnegative integer"))
    end
    top=isnothing(max_cell_dimension) ? full : min(full,max_cell_dimension)
    plan=certify_level_sections(M,certificate;faces,progress,on_section_certificate)
    parent=top==full ? full+1 : section_parent_dimension(plan,full+1,top)
    assemble_section_level(M,certificate,plan;faces,top,parent_max_dimension=parent,progress)
end

function assemble_section_level(M,certificate,plan;faces,top,parent_max_dimension,progress=false,extra_refinements=Tuple{Int,Int,Int}[])
    @assert certificate.passed
    C=M.complex;P=C.enumeration.problem;I=identity_matrix(P.field,P.rank)
    @assert 0<=top<parent_max_dimension<=length(C.cells)-1
    @assert all(r.parent[1]<=parent_max_dimension || r.parent[1]-r.core_dimension-2>=top for r in plan.proofs)
    sections=plan.sections;proofs=plan.proofs;lookup=plan.lookup
    marked=Set(t for t in union(plan.refined,Set(extra_refinements)) if t[1]<=parent_max_dimension+1)
    ids=[[sum(length,C.cells[1:k-1];init=0)+i for i in eachindex(layer)] for (k,layer) in enumerate(C.cells)]
    height(v)=certificate.bases[ids[v.dim+1][v.orbit]]+v.shift
    # Refining an entire original-cell section requires refining its boundary.
    # Conservatively include every crossing component in each original facet.
    for k in parent_max_dimension+1:-1:3
        for (layer,i,q) in collect(marked)
            layer==k || continue
            for (a,g) in C.boundaries[k][i]
                qD=q+Int(AV.unit_coordinates(P,M.units,det(g)).exponent)
                haskey(lookup[k-1][abs(a)],qD) && push!(marked,(k-1,abs(a),qD))
            end
        end
    end
    orbits=FlagOrbit[];flag_exponents=Int[]
    flag_lookup=Dict{Tuple{Int,Int,Int},Dict{Any,Tuple{Int,Any}}}()
    flag_groups=Dict{Tuple{Int,Int,Int},Vector{Int}}()
    spans=[[(minimum([height(M.vertices[k][i]);height.(faces[k][i])]),
             maximum([height(M.vertices[k][i]);height.(faces[k][i])])) for i in eachindex(layer)] for (k,layer) in enumerate(C.cells)]
    for (k,i,q) in sort!(collect(marked))
        table=get!(Dict{Any,Tuple{Int,Any}},flag_lookup,(k,i,q))
        group=get!(Vector{Int},flag_groups,(k,i,q))
        function descend(flag,S)
            h=height.(flag).+q
            if minimum(h)<0<maximum(h)
                id=length(orbits)+1
                if length(flag)<=top+1
                    for s in C.stabilizers[k][i]
                        key=Tuple(morse_act(M,s,v).key for v in flag)
                        @assert !haskey(table,key) || first(table[key])==id
                        table[key]=(id,s)
                    end
                end
                push!(orbits,FlagOrbit(flag,S));push!(flag_exponents,q);push!(group,id)
            end
            length(flag)>=top+2 && return
            v=first(flag);v.dim==0 && return
            lo,hi=spans[v.dim+1][v.orbit]
            min(minimum(h),lo+v.shift+q)<0<max(maximum(h),hi+v.shift+q) || return
            candidates=[morse_act(M,v.translate,w) for w in faces[v.dim+1][v.orbit]]
            seen=Set{Any}();scalar=all(s->s==I || s==-I,S)
            for w in candidates
                w.key in seen && continue
                next=scalar ? S : Any[]
                if !scalar
                    for s in S
                        image=morse_act(M,s,w);push!(seen,image.key)
                        image.key==w.key && push!(next,s)
                    end
                end
                descend([w;flag],next)
            end
        end
        descend([M.vertices[k][i]],C.stabilizers[k][i])
        clear_action_cache!(M)
    end
    cells=[SectionLevelCell[] for _ in 0:top]
    stabilizers=[Vector{Any}[] for _ in 0:top]
    cell_ids=Dict{Tuple{Symbol,Int,Int},Int}()
    translations=Dict{Tuple{Int,Int},Any}()
    flips=Dict{Tuple{Symbol,Int},Bool}()
    function add_cells(kind,id,d,q,S)
        flip=any(s->det(s)==-1,S);flips[(kind,id)]=flip
        for sgn in (flip ? (1,) : (1,-1))
            t=get!(translations,(q,sgn)) do
                diagonal_matrix(P.field,[sgn*M.units.epsilon^q;fill(one(P.field),P.rank-1)])
            end
            push!(cells[d+1],SectionLevelCell(kind,id,t,q,sgn))
            push!(stabilizers[d+1],[t*s*inv(t) for s in S if det(s)==1])
            cell_ids[(kind,id,sgn)]=length(cells[d+1])
        end
    end
    for (id,r) in enumerate(sections)
        (r.parent_dimension+1,r.parent_orbit,r.exponent) in marked && continue
        r.parent_dimension<=top+1 || continue
        add_cells(:section,id,r.parent_dimension-1,r.exponent,r.stabilizer)
    end
    for (id,r) in enumerate(orbits)
        add_cells(:flag,id,length(r.flag)-2,flag_exponents[id],r.stabilizer)
    end
    if progress
        println("Section level cells: ",length.(cells),"; refined original sections ",length(marked),"; flag orbits ",length(orbits));flush(stdout)
    end
    boundaries=[[Tuple{Int,Any}[] for _ in layer] for layer in cells]
    function transport(kind,id,k,t)
        S=kind==:section ? sections[id].stabilizer : orbits[id].stabilizer
        unit=AV.unit_coordinates(P,M.units,det(t))
        q=kind==:section ? sections[id].exponent : flag_exponents[id]
        @assert Int(unit.exponent)==q
        sgn=flips[(kind,id)] ? 1 : (-1)^unit.sign
        number=cell_ids[(kind,id,sgn)];target=cells[k][number]
        x=first(x for x in S if det(t)*det(x)==det(target.translate))
        transporter=t*x*inv(target.translate)
        @assert det(transporter)==1
        number,transporter,x
    end
    for k in 2:length(cells),(i,c) in enumerate(cells[k])
        if c.kind==:flag
            F=orbits[c.orbit].flag
            for j in eachindex(F)
                sub=[F[l] for l in eachindex(F) if l!=j]
                h=height.(sub).+c.exponent
                minimum(h)<0<maximum(h) || continue
                v=last(sub);q=c.exponent+v.shift
                key=Tuple(morse_act(M,inv(v.translate),w).key for w in sub)
                orbit,s=flag_lookup[(v.dim+1,v.orbit,q)][key]
                number,g,_=transport(:flag,orbit,k-1,c.translate*v.translate*s)
                augmentation=length(sub)==2 ? Int(sign(height(sub[2])-height(sub[1]))) : 1
                push!(boundaries[k][i],((-1)^j*augmentation*number,g))
            end
            continue
        end
        r=sections[c.orbit];parent=lookup[k+1][r.parent_orbit][r.exponent]
        for (a,g) in C.boundaries[k+1][r.parent_orbit]
            q=r.exponent+Int(AV.unit_coordinates(P,M.units,det(g)).exponent)
            haskey(lookup[k][abs(a)],q) || continue
            if (k,abs(a),q) in marked
                for id in flag_groups[(k,abs(a),q)]
                    F=orbits[id].flag;length(F)==k || continue
                    seen=Set{Any}()
                    for s in C.stabilizers[k][abs(a)]
                        image=[morse_act(M,s,v) for v in F]
                        key=Tuple(v.key for v in image);key in seen && continue;push!(seen,key)
                        transported=[morse_act(M,g,v) for v in image]
                        v=first(v for v in transported if parent.direction*(height(v)+r.exponent)>0)
                        parent.component[parent.positions[v.key]]==r.component || continue
                        number,t,_=transport(:flag,id,k-1,c.translate*g*s)
                        orientation=maximal_flag_orientation(M,image)
                        augmentation=k==2 ? Int(sign(height(image[2])-height(image[1]))) : 1
                        push!(boundaries[k][i],(-sign(a)*orientation*augmentation*number,t))
                    end
                end
            else
                lower=lookup[k][abs(a)][q];lower_top=M.vertices[k][abs(a)]
                for j in eachindex(lower.witnesses)
                    v=parent.direction*(height(lower_top)+q)>0 ? lower_top : lower.witnesses[j]
                    image=morse_act(M,g,v)
                    parent.component[parent.positions[image.key]]==r.component || continue
                    id,s=lower.maps[j]
                    number,t,x=transport(:section,id,k-1,c.translate*g*s)
                    orientation=AV.cell_orientation(P,C.cells[k][abs(a)],s*x)
                    augmentation=k==2 ? sections[id].augmentation : 1
                    push!(boundaries[k][i],(-sign(a)*orientation*augmentation*number,t))
                end
            end
        end
    end
    MorseSectionLevelComplex(M,certificate,proofs,sections,orbits,cells,stabilizers,boundaries,lookup,marked,parent_max_dimension)
end

function level_orientation(L::MorseSectionLevelComplex,k,i,g)
    c=L.cells[k][i];c.kind==:flag && return 1
    P=L.experiment.complex.enumeration.problem;r=L.sections[c.orbit]
    AV.cell_orientation(P,L.experiment.complex.cells[k+1][r.parent_orbit],inv(c.translate)*g*c.translate)
end


level_model(L::MorseSectionLevelComplex)=:sections
retained_parent_dimension(L::MorseSectionLevelComplex)=L.parent_dimension
