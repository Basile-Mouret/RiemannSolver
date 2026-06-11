using CairoMakie, GeometryBasics

"""
Plot 2D cell values
"""
function plot_cell_values(mesh::Mesh2D, cell_values::Vector{Float64}; title::String = "Finite Volume Plot")
    fig = Figure(fontsize=18, size=(600, 600))
    ax = Axis(fig[1, 1], title = title, aspect = DataAspect(), xlabel = "x", ylabel = "y")

    # 1. Build a discrete Polygon for each cell
    # This explicitly maps the exact coordinates to a closed shape, 
    # preventing Makie from scrambling or sharing vertices.
    polys = [
            Polygon([Point2f(mesh.points[v]) for v in c])
        for c in mesh.cells
    ]
    # 2. Plot the polygons
    # Setting strokewidth = 0 removes the triangle borders for a smooth heatmap look
    m = poly!(ax, polys, color = cell_values, colormap = :viridis, strokewidth = 1)

    # Add Colorbar
    Colorbar(fig[1, 2], m, label="Cell Values")
    
    return fig
end

"""
Animates a time dependent problem on a 2D unstructured mesh
"""
function animate_cell_values(
    mesh::Mesh2D,
    U_hist::Vector{Vector{Float64}},
    filename::String;
    dt_hist = nothing
)
    if !isnothing(dt_hist) && length(dt_hist) == length(U_hist)
        t_hist = cumsum(dt_hist)
    else
        t_hist = nothing
    end

    polys = [Polygon([Point2f(mesh.points[v]) for v in c]) for c in mesh.cells]

    cmin = minimum(minimum(u) for u in U_hist)
    cmax = maximum(maximum(u) for u in U_hist)
    crange = cmin == cmax ? (cmin - 1.0, cmax + 1.0) : (cmin, cmax)

    fig_anim = Figure(size = (700, 600))
    ax_anim = Axis(fig_anim[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
    obs_numerical = Observable(U_hist[1])
    m = poly!(ax_anim, polys, color = obs_numerical, colormap = :viridis, colorrange = crange, strokewidth = 1)
    Colorbar(fig_anim[1, 2], m)

    record(fig_anim, filename, 1:length(U_hist), framerate = 60) do frame_idx
        obs_numerical[] = U_hist[frame_idx]
        if !isnothing(dt_hist) && !isnothing(t_hist)
            t_val = t_hist[frame_idx]
            dt_val = dt_hist[frame_idx]
            ax_anim.title = "step $(frame_idx)  |  dt = $(round(dt_val, sigdigits=3))  |  t = $(round(t_val, digits=4))"
        else
            ax_anim.title = "Time Step: $(frame_idx - 1)"
        end
    end
    println("Saved animation to $filename")
end
