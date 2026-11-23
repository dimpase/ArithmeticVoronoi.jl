# Test the boundary poset on actual translates. Within a single closure, give
# each lifted cell (orbit, transporter coset) a local integer label.
function incidence_cycle(edges)
    isempty(edges) && return false
    degree = Dict{Int,Int}()
    for (u,v) in edges
        u == v && return false
        degree[u] = get(degree,u,0) + 1
        degree[v] = get(degree,v,0) + 1
    end
    all(==(2), values(degree)) || return false
    components, _ = section_components(keys(degree), edges)
    length(components) == 1
end

function regular_cell_boundary(boundary, canonical, k, i, identity)
    # k indexes the layer of cells of dimension k-1.
    2 <= k <= 4 || throw(ArgumentError("expected a cell of dimension one, two or three"))
    facets = Vector{Vector{Vector{Int}}}(undef,k)
    cells = Tuple{Int,Any}[(i,identity)]
    for layer in k:-1:2
        facets[layer] = Vector{Int}[]
        lower = Tuple{Int,Any}[]
        labels = Dict{Any,Int}()
        for (orbit,h) in cells
            terms = boundary(layer,orbit)
            isempty(terms) && return false
            if layer == 2
                length(terms) == 2 && sign(terms[1][1]) == -sign(terms[2][1]) || return false
            end
            ids = Int[]
            for (j,g) in terms
                key, rep, _ = canonical(layer-1,abs(j),h*g)
                id = get!(labels,(abs(j),key)) do
                    push!(lower,(abs(j),rep))
                    length(lower)
                end
                push!(ids,id)
            end
            allunique(ids) || return false
            push!(facets[layer],ids)
        end
        cells = lower
    end
    k == 2 && return true

    endpoints = facets[2]
    polygons = facets[3]
    all(edges -> incidence_cycle(endpoints[edges]), polygons) || return false
    k == 3 && return true

    # A connected closed surface with circular vertex links is a sphere
    # exactly when its Euler characteristic is two. Digonal faces are valid.
    length(cells) - length(endpoints) + length(polygons) == 2 || return false
    faces_at_edge = [Int[] for _ in endpoints]
    for (f,edges) in enumerate(polygons), e in edges
        push!(faces_at_edge[e],f)
    end
    all(fs -> length(fs) == 2, faces_at_edge) || return false
    components, _ = section_components(eachindex(polygons), faces_at_edge)
    length(components) == 1 || return false
    links = [Tuple{Int,Int}[] for _ in cells]
    for (e,vertices) in enumerate(endpoints)
        f1, f2 = faces_at_edge[e]
        for v in vertices
            push!(links[v],(f1,f2))
        end
    end
    all(incidence_cycle,links)
end

function preserves_regular_boundaries(proposed, boundaries, incident, canonical, identity)
    # Every new incidence occurs in a proposed row. All other parent
    # relations can be read from the current incidence index.
    to_check = collect(keys(proposed))
    seen = Set(to_check)
    for (k,i) in to_check
        # Higher-dimensional cells retain their certified geometric closures.
        k > 4 && return false
        for parent in incident[k][i]
            cell = (k+1,parent)
            if !(cell in seen)
                push!(seen,cell)
                push!(to_check,cell)
            end
        end
    end
    boundary(k,i) = get(proposed,(k,i),boundaries[k][i])
    all(cell -> regular_cell_boundary(boundary,canonical,cell...,identity), to_check)
end
