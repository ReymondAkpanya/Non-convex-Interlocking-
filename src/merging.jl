using LinearAlgebra
include("Polyhedron.jl")
include("affine_geometry")

"""
    merge!(poly1::Polyhedron, poly2::Polyhedron, facets1::Vector{<:Vector{<:Int}}, facets2::Vector{<:Vector{<:Int}})

Merge the two polyhedra poly1 and poly2 such that the corresponding vertices of the faces in faces1 and faces2 align. 
The aligned faces are deleted in the resulting polyhedron so that the result itself is a 2-dimensional piecewise linear manifold.
"""
function merge!(poly1::Polyhedron, poly2::Polyhedron, facets1::Vector{<:Vector{<:Int}}, facets2::Vector{<:Vector{<:Int}})
    @assert all([Set(f) in Set.(get_facets(poly1)) for f in facets1]) "faces1 needs to consist of facets of poly1."
    @assert all([Set(f) in Set.(get_facets(poly2)) for f in facets2]) "faces2 needs to consist of facets of poly2."

    @assert dimension(poly1) == dimension(poly2) "poly1 and poly2 need to be embedded into the same space."

    ind_to_align = sort(unique(hcat(vcat(facets1...), vcat(facets2...)), dims = 1), dims = 1)
    # display(ind_to_align)
    @assert length(unique(ind_to_align[:,1])) == length(ind_to_align[:,1]) "Merging isn't well defined"
    @assert length(unique(ind_to_align[:,1])) == length(ind_to_align[:,2]) "Merging isn't well defined"

    verts_to_align1 = get_verts(poly1)[ind_to_align[:,1]]
    verts_to_align2 = get_verts(poly2)[ind_to_align[:,2]]

    for i in 1:length(verts_to_align1)
        for j in i+1:length(verts_to_align2)
            @assert dist(verts_to_align1[i], verts_to_align1[j]) == dist(verts_to_align2[i], verts_to_align2[j]) "Polyhedra cannot be merged. Distance between Vertex $(ind_to_align[i,1]) and $(ind_to_align[j,1]) of poly1 is $(dist(verts_to_align1[i], verts_to_align1[j])), but the distance between vertex $(ind_to_align[i,2]) and $(ind_to_align[j,2]) of poly2 is $(dist(verts_to_align2[i], verts_to_align2[j]))"
        end
    end
    preim = verts_to_align2[affinebasis_indices(verts_to_align2)]
    im = verts_to_align1[affinebasis_indices(verts_to_align2)]
    @assert length(preim) >= 3 "Polyhedra cannot be merged along degenerate faces."

    if length(preim) == 3
        append!(preim, [preim[1] + cross(preim[2] - preim[1], preim[3] - preim[1])])
        append!(im, [im[1] - cross(im[2] - im[1], im[3] - im[1])])
    end

    tau = rigidmap(preim, im) # rigid map mapping facets2 to facets1
    
    # function that maps the index of a vertex of poly2 to the resulting vertex after merging
    function index_map2(i::Integer)
        if i in ind_to_align[:,2]
            return ind_to_align[indexin(i, ind_to_align[:,2])[1], 1]
        end

        sol = i + length(get_verts(poly1)) - sum(map(j -> i > j ? 1 : 0, ind_to_align[:,2]))
        return sol
    end

    sol_verts1 = get_verts(poly1)
    sol_verts2 = tau.(get_verts(poly2))[map(index -> !(index in ind_to_align[:,2]), 1:length(get_verts(poly2)))]
    sol_verts = vcat(sol_verts1, sol_verts2)

    sol_edges1 = get_edges(poly1)
    sol_edges2 = map(e -> index_map2.(e), get_edges(poly2))
    sol_edges = vcat(sol_edges1, sol_edges2)

    sol_facets1 = get_facets(poly1)[map(f -> !(Set(f) in Set.(facets1)), get_facets(poly1))]
    sol_facets2 = get_facets(poly2)[map(f -> !(Set(f) in Set.(facets1)), get_facets(poly2))]
    sol_facets2 = [index_map2.(f) for f in sol_facets2]
    sol_facets = vcat(sol_facets1, sol_facets2)

    set_verts!(poly1, sol_verts)
    set_edges!(poly1, sol_edges)
    set_facets!(poly1, sol_facets)
end


"""
    merge(poly1::Polyhedron, poly2::Polyhedron, facets1::Vector{<:Vector{<:Int}}, facets2::Vector{<:Vector{<:Int}})

Merge the two polyhedra poly1 and poly2 such that the corresponding vertices of the faces in faces1 and faces2 align. 
The aligned faces are deleted in the resulting polyhedron so that the result itself is a 2-dimensional piecewise linear manifold.
"""
function merge(poly1::Polyhedron, poly2::Polyhedron, facets1::Vector{<:Vector{<:Int}}, facets2::Vector{<:Vector{<:Int}})
    poly1_copy = deepcopy(poly1)
    merge!(poly1_copy, poly2, facets1, facets2)

    return poly1_copy
end