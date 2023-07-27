using PlotlyJS
using Colors
import PlotlyJS.plot

include("Polyhedron.jl")
include("decomposition.jl")

"""
Aux function to plot a polyhedron. Returns array of traces that can be handled by PlotlyJS.
"""
function trace_polyhedron(poly::Polyhedron; color::Color = RGB(0,0.9,1), labels::Bool = false, opacity::Real = 0.6)
    polyTriang = triangulate(poly)

    facecolor = repeat([color], length(polyTriang.facets))

    # plot faces
    mesh = mesh3d(
        x = [vert[1] for vert in polyTriang.verts],
        y = [vert[2] for vert in polyTriang.verts],
        z = [vert[3] for vert in polyTriang.verts],
        i = [triang[1] for triang in polyTriang.facets].-1,
        j = [triang[2] for triang in polyTriang.facets].-1,
        k = [triang[3] for triang in polyTriang.facets].-1,
        facecolor = facecolor,
        opacity = opacity
    )

    traces = [mesh]

    # plot edges
    for edge in poly.edges
        trace = scatter(
            x = [poly.verts[v][1] for v in edge],
            y = [poly.verts[v][2] for v in edge],
            z = [poly.verts[v][3] for v in edge],
            line = attr(color = "black"),
            mode = "lines",
            type = "scatter3d"
        )

        push!(traces, trace)
    end

    # plot vertices
    mode = labels ? "markers+text" : "markers"
    trace = scatter3d(
    	x = [v[1] for v in poly.verts],
        y = [v[2] for v in poly.verts],
        z = [v[3] for v in poly.verts],
        mode = mode,
        text = [string(v) for v in 1:length(poly.verts)],
        textposition = "bottom center",
        marker = attr(color = color, size = 4)
    )

    push!(traces, trace)

    return traces
end


function plot(assembly::Vector{Polyhedron}; labels::Bool = false, width::Int = 600, height::Int = 600)
    colors = distinguishable_colors(length(assembly), RGB(0,0,0))

    plot(reverse(union([trace_polyhedron(poly; color = colors[i], labels = labels) for (i, poly) in enumerate(assembly)]...)), 
         Layout(showlegend = false, 
                autosize = false, 
                width = width, 
                height = height,
                scene_aspectmode = "data"
         ))
end


function plot(poly::Polyhedron; color::Color = RGB(0,0.9,1), labels::Bool = false, width::Int = 600, height::Int = 600)
    plot(trace_polyhedron(poly; color = color, labels = labels), 
         Layout(showlegend = false, autosize = false, width=width, height = height, scene_aspectmode = "data"))
end