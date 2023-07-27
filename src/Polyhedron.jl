using GenericLinearAlgebra
using PolygonOps

include("affine_geometry")

abstract type AbstractPolyhedron end
mutable struct Polyhedron <:AbstractPolyhedron
    verts::Vector{<:Vector{<:Real}} # vertex array. Every vertex is an array of 3 spatial coordinates
    edges::Vector{<:Vector{<:Integer}} # edge array. Every edge is an array of the indices of the adjacent vertices
    facets::Vector{<:Vector{<:Integer}} # facet array. Every facet is an array of the indices on its boundary. The last vertex is adjacent to the first.
    # TODO: Neuen constructor mit optionalen Argumenten (wenn nur coordinates gegeben werden, ist Ergebnis die konvexe Hülle der Punkte + check der Dimension
    # wenn nur Facets gegeben sind, werden Edges automatisch gesetzt und es wird gecheckt, dass Vertizes auf einer Facet koplanar aber nicht kollinear sind)
end

function get_verts(poly::Polyhedron)
    return deepcopy(poly.verts)
end

function set_verts!(poly::Polyhedron, verts::Vector{<:Vector{<:Real}})
    d = length(verts[1])
    @assert d == 3 "Only 3-dimensional polyhedra supported."
    @assert all([length(v) == d for v in verts]) "Dimension mismatch in vertices."
    poly.verts = verts
end

function get_edges(poly::Polyhedron)
    return deepcopy(poly.edges)
end

function set_edges!(poly::Polyhedron, edges::Vector{<:Vector{<:Int}})
    @assert all(length(e) == 2 for e in edges) "Edges need to consist of vectors of length 2."
    # @assert sort(union(edges...)) == [1:max(union(edges...)...)...] "Vertex indices need to be 1, ..., $(length(unique(vcat(edges...))))."
    # TODO: Assert, dass die Kanten auch tatsächlich auf Rand von Facets liegen?
    poly.edges = edges
end


function get_facets(poly::Polyhedron)
    return deepcopy(poly.facets)
end

function set_facets!(poly::Polyhedron, facets::Vector{<:Vector{<:Int}})
    # @assert sort(union(facets...)) == [1:max(union(facets...)...)...] "Vertex indices need to be 1, ..., $(length(unique(vcat(facets...))))."
    @assert all([affinedim(get_verts(poly)[f]) == 2 for f in facets]) "Facets have to span affine spaces of dimension 2."
    poly.facets = facets
end

"""
    dimension(poly::Polyhedron)

Get the dimension of the unerlying space the polyhedron is embedded into.
"""
function dimension(poly::Polyhedron)
    return length(get_verts(poly)[1])
end


"""
    isadjacent(poly::Polyhedron, facetoredge::Vector{<:Int}, facet::Vector{<:Int})


Calculate whether the facet or edge facetoredge and the facet facet of the polyhedron poly are adjacent, i.e. share a common edge.
"""
function isadjacent(poly::Polyhedron, facetoredge::Vector{<:Int}, facet::Vector{<:Int})
    @assert Set(facetoredge) in Set.(get_facets(poly)) || Set(facetoredge) in Set.(get_edges(poly)) "facetoredge has to be a facet or an edge of poly."
    @assert Set(facet) in Set.(get_facets(poly)) "facet2 has to be a facet of poly."

    intersection = Base.intersect(facetoredge, facet)

    return any(map(edge -> Base.intersect(edge, intersection) == edge, get_edges(poly)))
end


"""
    adjfacets(poly::Polyhedron, facetoredge::Vector{<:Int})

Calculate the adjacent facets of the facet or edge facet in the Polyhedron poly, i.e. the facets sharing at least one edge with facet.
"""
function adjfacets(poly::Polyhedron, facetoredge::Vector{<:Int})
    return setdiff(get_facets(poly)[map(facet2 -> isadjacent(poly, facetoredge, facet2), get_facets(poly))], [facetoredge])
end

"""
    formpath!(vertexindices::Vector{<:Int}, poly::Polyhedron)

Sort the vector vertexindices such that the corresponding vertices of the Polyhedron poly form a vertex edges path.
"""
function formpath!(vertexindices::Vector{<:Int}, poly::Polyhedron)
    endpoints = vertexindices[map( i -> length( Base.intersect([Set([i,j]) for j in vertexindices], Set.(get_edges(poly))) ) == 1, vertexindices )]
    @assert length(endpoints) <= 2 "Vertices don't lie on a common path"
    intersectionverts = vertexindices[map( i -> length( Base.intersect([Set([i,j]) for j in vertexindices], Set.(get_edges(poly))) ) > 2, vertexindices )]
    @assert length(intersectionverts) == 0 "No intersections allowed."

    tosort = deepcopy(vertexindices)
    if length(endpoints) == 0
        start = vertexindices[1]
    else
        start = endpoints[1]
    end
    path = [start]
    setdiff!(tosort, [start])
    vertexindices[1] = start

    for i =2:length(vertexindices)
        next = tosort[map( j -> (Set([vertexindices[i-1], j]) in Set.(get_edges(poly))), tosort )][1]
        vertexindices[i] = next
        setdiff!(tosort, [next])
    end
end

function formpath(vertexindices::Vector{<:Int}, poly::Polyhedron)
    vertexindicescopy = deepcopy(vertexindices)
    formpath!(vertexindicescopy, poly)
    return vertexindicescopy
end

struct Ray
    point::Vector{<:Real}
    vector::Vector{<:Real} 
end

struct Plane
    point::Vector{<:Real} 
    vectors::Vector{<:Vector{<:Real}} 
end


"""
A::Matrix
cond::Float64
Returns an orthonormal basis for the columnspace of the matrix A using svd. Singular values with abs < cond are treated as 0.
"""
function colspace(A::Matrix{<:Real}; tol::Real = 1e-5)
    F = svd(A)
    return [c[:] for c in eachcol(F.U[:, findall(>(tol), abs.(F.S))])]
end


"""
A::Matrix
cond::Float64
Returns an orthonormal basis for the nullspace of the matrix A using svd. Singular values with abs < cond are treated as 0.
"""
function nullspace(A::Matrix{<:Real}, tol::Real = 1e-5)
    F = svd(A)
    return [c[:] for c in eachcol(F.V[:, findall(<(tol), abs.(F.S))])]
end

"""
returns the plane, in which polygon lies.
"""
function plane(polygon::Vector{<:Vector{<:Real}}; tol::Real=1e-8)::Plane 
    @assert polygon[1] == polygon[end] "First and last vertex of polygon have to equal."

    basis = affinebasis(polygon, atol = tol)
    @assert length(basis) == 3 "polygon doesn't span a plane."
    v = basis[1]
    A = map(b -> b-v, basis[2:end])
    return Plane(v, A)
end

"""
returns the intersection between a ray and a plane if it exists.
"""
function intersect(ray::Ray, plane::Plane; tol::Real = 1e-8)
    vr = ray.point
    vp = plane.point

    A = hcat(union(plane.vectors, [-ray.vector])...)
    if length(colspace(A, tol=tol)) < 3
        throw(ErrorException("ray and plane are parallel."))
    end
    
    coeffs = A\(vr - vp)
    if coeffs[3] < 0
        throw(ErrorException("ray and plane don't intersect"))
    end

    return vr + coeffs[3] * ray.vector
end


"""
returns 1 if the 3d point p lies in the polygon poly embedded into R^3, -1 if it lies on its boundary and 0 otherwise
"""
function inpolygon3d(p::AbstractVector, poly::AbstractVector; tol::Real=1e-5)
    d = 3
    @assert length(p) == 3 "Only 3d case."
    @assert all([length(v) == 3 for v in poly]) "Only 3d case."

    E = plane(poly, tol = tol)
    point = deepcopy(p)

    preim = [E.point, E.point + E.vectors[1], E.point + E.vectors[2], E.point + normalize(cross(E.vectors[1], E.vectors[2]))]
    im = [[0,0,0], [1,0,0], [0,1,0], [0,0,1]]

    # affine map so that the affine basis of the polygon is mapped to the xy plane
    aff = affinemap(preim, im, atol = tol)

    # transformed polygon and point
    polytransformed = aff.(poly)
    pointtransformed = aff(point)
    

    if abs(pointtransformed[3]) < tol
         # inpolygon projects onto the xy plane and decides the problem.
         # If p and poly lie in the same plane, then this is enough to decide the problem
        return inpolygon(pointtransformed, polytransformed)
    else
        # p and poly don't lie in the same plane
        return 0 
    end
end


# decide whether point is inside of polyhedron
"""
    Randomized algorithm to check whether a point is contained in a polyhedron.
"""
function inpolyhedron(point::Vector{<:Real}, poly::Polyhedron; tol::Real=1e-5)::Int 
    # check whether point lies on the boundary of poly
    for facet in get_facets(poly)
        polygon = push!(map(v -> poly.verts[v], facet), poly.verts[facet[1]])
        if  inpolygon3d(point, polygon, tol=tol) != 0
            return -1
        end
    end

    while true
        r = Ray(point, normalize!(randn(Float64, 3)))
        numIntersections = 0 # number of intersections of r and the facets of poly

        for facet in get_facets(poly)
            polygon = push!(map(v -> poly.verts[v], facet), poly.verts[facet[1]])
            E = plane(polygon)

            try
                intersect(r, E)
            catch error
                continue
            end

            p = intersect(r, E)

            if inpolygon3d(p, polygon, tol = tol) == -1
                break
            elseif inpolygon3d(p, polygon, tol = tol) == 1
                numIntersections = numIntersections + 1
            end
        end

        if mod(numIntersections, 2) == 0
            return 0
        else
            return 1
        end
    end
end



