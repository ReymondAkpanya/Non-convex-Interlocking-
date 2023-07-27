include("Polyhedron.jl")
include("combinatorics.jl")


# triangulate surface of Polyhedron
"""
returns a polyhedron containing the vertices and edges of poly such that every facet is triangular.
"""
function triangulate(poly::Polyhedron)::Polyhedron
    newVerts = deepcopy(poly.verts)
    newEdges = deepcopy(poly.edges)
    newFacets = Vector{Int}[]


    for facet in poly.facets
        subfacet = deepcopy(facet)
        polygon = push!(map(i -> poly.verts[i], subfacet), poly.verts[subfacet[1]])

        while length(subfacet) > 3
            if inpolygon3d((poly.verts[subfacet[end]] + poly.verts[subfacet[2]]) / 2, polygon) == 1
                push!(newEdges, [subfacet[end], subfacet[2]])
                push!(newFacets, [subfacet[end], subfacet[1], subfacet[2]])
                subfacet = subfacet[2:end]
                continue
            elseif inpolygon3d((poly.verts[subfacet[end - 1]] + poly.verts[subfacet[1]]) / 2, polygon) == 1
                push!(newEdges, [subfacet[end - 1], subfacet[1]])
                push!(newFacets, [subfacet[end - 1], subfacet[length(subfacet)], subfacet[1]])
                subfacet = subfacet[1:end - 1]
                continue
            else
                for ind in 2:length(subfacet) - 1
                    if inpolygon3d((poly.verts[subfacet[ind - 1]] + poly.verts[subfacet[ind + 1]]) / 2, polygon) == 1
                        push!(newEdges, [subfacet[ind - 1], subfacet[ind + 1]])
                        push!(newFacets, [subfacet[ind - 1], subfacet[ind], subfacet[ind + 1]])
                        subfacet = subfacet[1:end .!= ind]
                        break
                    end
                end
            end
        end

        push!(newFacets, subfacet)
    end

    return Polyhedron(newVerts, newEdges, newFacets)
end

"""
e::Vector{Int}
poly::Polyhedron
tol::Float64
Checks whether the edge e is turnable edge of poly. I.e. poly is 
    - a spherical simplicial surface and 
    - the line connecting the wingtips of the butterfly with inner edge e is contained in poly.
Floats with abs value <tol are considered zeros
"""
function isturnable(e::Vector{<:Int}, polyhedron::Polyhedron, tol::Real=1e-5)::Bool
    @assert all(length.(polyhedron.facets).==3) "poly may only have triangle facets"
    @assert in(e, polyhedron.edges) || in(reverse(e), poly.edges) "e has to be an edge of poly"
    
    butterfly = polyhedron.facets[map(f -> length(Base.intersect(f, e)) == 2, polyhedron.facets)]
    butterflyTips = setdiff(union(butterfly...), e)
    if inpolyhedron((polyhedron.verts[butterflyTips[1]] + polyhedron.verts[butterflyTips[2]]) / 2, polyhedron, tol) == 0
        return 0
    end

    return 1
end


"""
poly::Polyhedron
tol::Float64
returns whether poly is convex
Floats with abs value <tol are considered zeros
"""
function isconvex(poly::Polyhedron, tol::Real=1e-5)::Bool
    polyhedron = deepcopy(poly)
    if any(length.(polyhedron.facets).!=3)
        polyhedron = triangulate(polyhedron)
    end

    for edge in polyhedron.edges
        if !isturnable(edge, polyhedron, tol)
            return 0
        end
    end
    return 1
end


"""
poly::Polyhedron
returns a vector of tetrahedra, which union is the Polyhedron poly.
"""
function convexdecomp(poly::Polyhedron)::Vector{Polyhedron}
    sol = Polyhedron[]

    triangPoly = triangulate(poly)
    subPoly = deepcopy(triangPoly)

    while !isconvex(subPoly)
        vertexDegrees = map(v -> length(FacesOfVertex(subPoly, v)), 1:length(subPoly.verts))        

        # remove a tetrahedron from subPoly if one exists
        if any(vertexDegrees.==3)
            # println("\nsubPoly:", subPoly)
            # println("vertexDegrees.==3:", [vertexDegrees.==3])
            v = indexin(3, vertexDegrees)[1]
            
            tetraVerts = push!(union(map(e -> setdiff(e, [v]), EdgesOfVertex(subPoly, v))...), v) # determine vertex indices of tetrahedron
            # println("tetraverts:", tetraVerts)
            tetra = Polyhedron(map(w -> subPoly.verts[w], tetraVerts), [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]], [[1,2,3], [1,2,4], [1,3,4], [2,3,4]]) # resulting tetrahedron
            # println("found tetra:", tetra)
            push!(sol, tetra)

            # remove v and adjacent edges and facets from subPoly
            subPoly.verts = subPoly.verts[1:end .!= v]
            setdiff!(subPoly.edges, subPoly.edges[map(e -> v in e, subPoly.edges)])
            setdiff!(subPoly.facets, subPoly.facets[map(f -> v in f, subPoly.facets)])

            # add bottom of tetrahedron as facet
            push!(subPoly.facets, tetraVerts[1:3])

            # shift indices of all vertices with index > v
            subPoly.edges = map(e -> map((w -> w>v ? w-1 : w), e), subPoly.edges)
            subPoly.facets = map(f -> map((w -> w>v ? w-1 : w), f), subPoly.facets)

            # println(sol)
            continue
        end

        # edge turn if no tetrahedron was removed
        for edge in subPoly.edges
            # if edge is turnable, cut the resulting tetrahedron out
            if isturnable(edge, subPoly)
                butterfly = subPoly.facets[map(f -> length(Base.intersect(f, edge)) == 2, subPoly.facets)]
                butterflyTips = setdiff(union(butterfly...), edge)

                # println("\nsubpoly: ", subPoly)
                # println("evaluated edge: ", edge)
                # println("corresponding butterfly: ", butterfly)
                setdiff!(subPoly.edges, [edge])
                push!(subPoly.edges, butterflyTips)

                setdiff!(subPoly.facets, butterfly)
                append!(subPoly.facets, [[edge[1], butterflyTips[1], butterflyTips[2]], [edge[2], butterflyTips[1], butterflyTips[2]]])

                tetra = Polyhedron(subPoly.verts[[edge[1], edge[2], butterflyTips[1], butterflyTips[2]]], [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]], [[1,2,3], [1,2,4], [1,3,4], [2,3,4]])
                push!(sol, tetra)
                # println("length(sol) = ", length(sol))
                break
            end
        end
    end

    push!(sol, subPoly) # subPoly is now convex

    return sol
end