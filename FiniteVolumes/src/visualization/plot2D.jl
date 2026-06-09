
function show_heatmap(mesh::Mesh2D, cellValues::Vector{Float64}, title::String)
    fig = Figure(fontsize=18, size=(600, 600))
    ax = Axis(fig[1, 1], title = title, aspect = DataAspect(), xlabel = "x", ylabel = "y")

    # 1. Build a discrete Polygon for each cell
    # This explicitly maps the exact coordinates to a closed shape, 
    # preventing Makie from scrambling or sharing vertices.
    polys = [
        Polygon([
            Point2f(mesh.points[c[1]]), 
            Point2f(mesh.points[c[2]]), 
            Point2f(mesh.points[c[3]])
        ]) 
        for c in mesh.cells
    ]
    # 2. Plot the polygons
    # Setting strokewidth = 0 removes the triangle borders for a smooth heatmap look
    m = poly!(ax, polys, color = cellValues, colormap = :viridis, strokewidth = 1)

    # Add Colorbar
    Colorbar(fig[1, 2], m, label="Cell Values")
    
    return fig
end

