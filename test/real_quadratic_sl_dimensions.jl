using ArithmeticVoronoi, Oscar, Test, JSON
const RQSL=ArithmeticVoronoi.RealQuadraticSL

@testset "dimension-independent determinant-height construction" begin
    # An eight-dimensional simplex goes beyond every rank-two link dimension.
    facets=[Tuple(1:9)]
    collapse=RQSL.collapse_link(facets)
    @test collapse.passed
    @test RQSL.replay_collapse(facets,collapse.pairs)
    @test_throws AssertionError RQSL.replay_collapse(facets,collapse.pairs[2:end])
    sphere=[Tuple(i for i in 1:6 if i!=j) for j in 1:6]
    @test !RQSL.collapse_link(sphere).passed
    @test isempty(RQSL.subset_masks(0))
    @test collect(RQSL.subset_masks(0;include_empty=true))==[0]
    @test collect(Iterators.take(RQSL.subset_masks(70),3))==BigInt[1,2,3]
    @test last(RQSL.subset_masks(70))==(big(1)<<70)-1
    @test_throws ArgumentError RQSL.subset_masks(-1)
    @test_throws AssertionError RQSL.certify_morse_links(Dict("base_heights"=>[0.0],"links"=>Any[]))
    for n in (2,3,4)
        P=arithmetic_problem(5;rank=n)
        units=ArithmeticVoronoi.verified_units(P)
        I=identity_matrix(P.field,n)
        one_step=diagonal_matrix(P.field,[units.epsilon;fill(one(P.field),n-1)])
        for (g,shift) in ((one_step,1),(units.epsilon*I,n))
            image=ArithmeticVoronoi.act_form(P,I,inv(g))
            @test RQSL.morse_base_height(P,units,image)≈RQSL.morse_base_height(P,units,I)+shift
        end
        # The full real SL centre has order two exactly in even rank.
        V=shortest_data(trace_gram(P,I)).vectors
        stab=ArithmeticVoronoi.cell_maps(P,V,V;determinant_one=true)
        @test length(stab)==factorial(n)*2^(n-1)
        @test all(g->det(g)==1,stab)
        @test (-I in stab)==iseven(n)
        @test ArithmeticVoronoi.hap_request(voronoi_request(5;group=:SL,rank=n))["rank"]==n
    end
end

@testset "flag orbits cover the independently enumerated face lattice" begin
    P=arithmetic_problem(5)
    C=voronoi_complex(enumerate_perfect_forms(P;backend=:native))
    M=RQSL.morse_experiment(C)
    raw=[[Vector{RQSL.MorseVertex}[] for c in layer] for layer in C.cells]
    key(flag)=Tuple(v.key for v in flag)
    for (k,layer) in enumerate(C.cells),(i,c) in enumerate(layer)
        top=M.vertices[k][i]
        if k==1
            raw[k][i]=[[top]]
        else
            for (j,g) in C.boundaries[k][i],flag in raw[k-1][abs(j)]
                push!(raw[k][i],[RQSL.MorseVertex[RQSL.morse_act(M,g,v) for v in flag];top])
            end
        end
        expanded=Set(key([RQSL.morse_act(M,s,v) for v in flag])
            for flag in M.orbit_flags[k][i] for s in C.stabilizers[k][i])
        @test expanded==Set(key.(raw[k][i]))
        @test length(expanded)==M.full_flag_counts[k][i]
        if k>1
            exhaustive=Set{Any}()
            for lower in C.cells[k-1],h in ArithmeticVoronoi.cell_maps(P,c.vectors,lower.vectors)
                push!(exhaustive,ArithmeticVoronoi.cell_vector_key(ArithmeticVoronoi.cell_image(P,inv(h),lower.vectors)))
            end
            actual=Set(ArithmeticVoronoi.cell_vector_key(ArithmeticVoronoi.cell_image(P,g,C.cells[k-1][abs(j)].vectors)) for (j,g) in C.boundaries[k][i])
            @test actual==exhaustive
        end
    end
end

@testset "exact graph-collapse witnesses" begin
    edges=[[i,j] for i in 1:9 for j in i+1:9]
    result=RQSL.collapse_flag_graph(1:9,edges)
    @test result.passed
    @test RQSL.replay_flag_collapse(1:9,edges,result.steps)
    @test_throws AssertionError RQSL.replay_flag_collapse(1:9,edges,result.steps[2:end])
    @test !RQSL.collapse_flag_graph(1:4,[[1,2],[2,3],[3,4],[4,1]]).passed
    @test_throws AssertionError RQSL.replay_flag_collapse(1:4,[[1,2],[2,3],[3,4],[4,1]],[[1,2]])
    triangle=[[1,2],[1,3],[2,3]]
    steps=[[1,2,3],[1,3],[2,3]]
    @test RQSL.replay_flag_collapse(1:3,triangle,steps)
    # Independently expand the edge-collapse macro into elementary face pairs.
    @test RQSL.replay_collapse([(1,2,3)],[((1,2),(1,2,3)),((1,),(1,3)),((2,),(2,3))])
    @test_throws AssertionError RQSL.replay_flag_collapse(1:3,[[1,2],[2,3]],steps)
end

@testset "incidence-poset links equal full flag links" begin
    P=arithmetic_problem(5)
    M=RQSL.morse_experiment(voronoi_complex(enumerate_perfect_forms(P;backend=:native)))
    legacy,oldlinks=RQSL.morse_link_data(M)
    data,newlinks,faces=RQSL.morse_link_graphs(M)
    @test data["base_heights"]==legacy["base_heights"]
    io=IOBuffer();JSON.json(io,data)
    saved=JSON.parse(String(take!(io)),Dict{String,Any})
    @test length(saved["links"])==length(oldlinks)
    @test length(saved["links"][1]["edges"])==length(data["links"][1]["edges"])
    for i in eachindex(oldlinks)
        old=oldlinks[i]
        new=sort!(collect(values(newlinks[i]));by=v->v.key)
        @test Set(v.key for v in old.vertices)==Set(v.key for v in new)
        edgekey(a,b)=minmax(a.key,b.key)
        expected=Set(edgekey(old.vertices[f[a]],old.vertices[f[b]]) for f in old.facets for a in eachindex(f) for b in a+1:length(f))
        actual=Set(edgekey(new[a],new[b]) for (a,b) in data["links"][i]["edges"])
        @test actual==expected
    end
    @test RQSL.certify_morse_graphs(data).passed
    # Interning must preserve an exact transporter and its determinant height.
    for layer in faces,face_list in layer,v in face_list
        actual=ArithmeticVoronoi.cell_image(P,v.translate,M.complex.cells[v.dim+1][v.orbit].vectors)
        @test RQSL.morse_vertex_key(v.dim,actual)==v.key
        @test v.shift==ArithmeticVoronoi.unit_coordinates(P,M.units,det(v.translate)).exponent
    end
    # Large unit powers must stay exact when compact machine integers no longer fit.
    g=diagonal_matrix(P.field,[M.units.epsilon^100,one(P.field)])
    v=M.vertices[1][1]
    image=RQSL.morse_act(M,g,v)
    @test image.shift==100
    @test image.key==RQSL.morse_vertex_key(v.dim,ArithmeticVoronoi.cell_image(P,g,M.complex.cells[1][1].vectors))
    @test any(abs(x)>typemax(Int) for ray in image.key[2] for x in ray)
    half=diagonal_matrix(P.field,[M.units.epsilon^50,one(P.field)])
    @test RQSL.morse_act(M,half,RQSL.morse_act(M,half,v)).key==image.key
    # Both factors fit machine integers, but their product requires BigInt.
    factor=diagonal_matrix(P.field,[M.units.epsilon^60,one(P.field)])
    first_image=RQSL.morse_act(M,factor,v)
    @test all(x isa Int for ray in first_image.key[2] for x in ray)
    product=RQSL.morse_act(M,factor,first_image)
    @test product.shift==120
    @test product.key==RQSL.morse_vertex_key(v.dim,ArithmeticVoronoi.cell_image(P,factor^2,M.complex.cells[1][1].vectors))
end

@testset "exact regular-level selection" begin
    # Sweep predictions agree with an independent direct count in every
    # interval between critical heights, including translated edge endpoints.
    bases=BigInt[1,3,7].//10
    edges=[(1,2,0,2),(2,3,1,1),(3,1,-1,2)]
    offset,before,after=RQSL.least_vertex_level(bases,edges)
    function count_at(t)
        sum(edges) do (i,j,q,w)
            lo,hi=minmax(bases[i],bases[j]+q)
            w*count(k->lo+k<t<hi+k,-5:5)
        end
    end
    candidates=[big(0)//1,1//5,1//2,17//20]
    @test before==count_at(0)
    @test after==count_at(offset)==minimum(count_at,candidates)
    @test all(!isinteger,bases .- offset)
    @test [a-b for a in bases for b in bases]==[a-b for a in bases.-offset for b in bases.-offset]
    P=arithmetic_problem(5)
    M=RQSL.morse_experiment(voronoi_complex(enumerate_perfect_forms(P;backend=:native));maximal_flags=false)
    data,_,faces=RQSL.morse_link_graphs(M)
    certificate=RQSL.certify_morse_graphs(data)
    L=RQSL.build_level_complex(M,data,certificate;max_cell_dimension=0,faces)
    @test length(only(L.cells))==data["level_selection"]["zero_cells_after"]
end

@testset "dimension-preserving symmetry collapses" begin
    edges=[(1,2),(2,3),(3,4),(4,1),(1,5),(2,5),(3,5),(4,5)]
    colours=[0,1,0,1,2]
    generators=[[3,4,1,2,5]]
    result=RQSL.orbit_collapse_flag_graph(1:5,edges,generators,colours)
    @test result.passed
    @test RQSL.replay_orbit_flag_collapse(1:5,edges,generators,colours,result.steps)
    @test result.expanded_steps>length(result.steps)
    saved=JSON.parse(JSON.json(result.steps))
    @test RQSL.replay_orbit_flag_collapse(1:5,edges,generators,colours,saved)
    @test_throws AssertionError RQSL.replay_orbit_flag_collapse(1:5,edges,generators,colours,result.steps[2:end])
    # Independently expand each orbit into individually checked graph steps.
    expanded=Tuple[]
    for step in result.steps
        queue=[step];seen=Set([step[1:end-1]])
        for entry in queue,g in generators
            next=map(v->g[v],entry);face=next[1:end-1]
            if !(face in seen);push!(seen,face);push!(queue,next);end
        end
        append!(expanded,queue)
    end
    @test RQSL.replay_flag_collapse(1:5,edges,expanded)
    # Proper colours alone are insufficient: this permutation is not a graph
    # automorphism, and the independent replay must reject it.
    @test_throws AssertionError RQSL.replay_orbit_flag_collapse(1:5,[(1,2),(2,3),(3,4),(4,5)],
        [[3,2,1,4,5]],[0,1,0,1,2],Tuple[])
    @test_throws ArgumentError RQSL.replay_morse_graph_certificate(Dict(),Int[],(model=:unknown,steps=Tuple[]))
    @test_throws AssertionError RQSL.orbit_collapse_flag_graph(1:5,edges,[[2,3,4,1,5]],colours)
    # An orbit of two cone edges expands to ordinary valid removals.
    triangles=[(1,2),(1,3),(2,3),(4,5),(4,3),(5,3)]
    @test RQSL.replay_flag_collapse(1:5,triangles,[(1,2,3),(4,5,3),(1,3),(2,3),(4,3),(5,3)])
end

@testset "collapsible links beyond cones" begin
    path=[(1,2),(2,3),(3,4)]
    edges=[path;[(v,a) for a in (5,6) for v in 1:4]]
    generators=[[1,2,3,4,6,5]]
    colours=[0,1,0,1,2,2]
    inner=[(1,2),(2,3),(3,4)]
    steps=Any[Dict("vertex"=>5,"link_steps"=>inner);inner]
    @test RQSL.replay_orbit_flag_collapse(1:6,edges,generators,colours,steps)
    @test RQSL.replay_orbit_flag_collapse(1:6,edges,generators,colours,JSON.parse(JSON.json(steps)))
    bad=Any[Dict("vertex"=>5,"link_steps"=>inner[2:end]);inner]
    @test_throws AssertionError RQSL.replay_orbit_flag_collapse(1:6,edges,generators,colours,bad)
    # Expand the relative star collapses independently into elementary pairs.
    pairs=Any[]
    for a in (5,6)
        append!(pairs,[((v,a),(v,v+1,a)) for v in 1:3])
        push!(pairs,((a,),(4,a)))
    end
    append!(pairs,[((v,),(v,v+1)) for v in 1:3])
    facets=[(v,v+1,a) for a in (5,6) for v in 1:3]
    @test RQSL.replay_collapse(facets,pairs)
    # A circular link cannot be used to delete the vertex of a cone.
    cycle=vcat(path,[(1,4)],[(v,5) for v in 1:4])
    @test_throws AssertionError RQSL.replay_orbit_flag_collapse(1:5,cycle,Vector{Int}[],colours[1:5],
        Any[Dict("vertex"=>5,"link_steps"=>inner)])
end

@testset "proper colours filter all possible cone apices" begin
    # Exhaust every graph with these proper colours. Compare the candidate
    # filter against the full definition of a vertex/edge cone, without using
    # colour counts in the independent test.
    colours=[0,1,2,0,1,2]
    possible=[(i,j) for i in 1:6 for j in i+1:6 if colours[i]!=colours[j]]
    for mask in 0:(1<<length(possible))-1
        edges=[e for (k,e) in enumerate(possible) if !iszero(mask & (1<<(k-1)))]
        graph=RQSL.flag_adjacency(1:6,edges)
        for v in 1:6
            full=Set(w for w in graph[v] if w!=v && issubset(graph[v],graph[w]))
            filtered=Set(w for w in RQSL.unique_colour_candidates(graph[v],colours)
                if w!=v && issubset(graph[v],graph[w]))
            @test full==filtered
        end
        for (u,v) in edges
            common=setdiff(intersect(graph[u],graph[v]),Set([u,v]))
            full=Set(w for w in common if issubset(common,graph[w]))
            filtered=Set(w for w in RQSL.unique_colour_candidates(common,colours) if issubset(common,graph[w]))
            @test full==filtered
        end
    end
end
