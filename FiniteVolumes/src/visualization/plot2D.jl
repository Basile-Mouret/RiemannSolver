using GLMakie, GeometryBasics
import CairoMakie
GLMakie.activate!()

"""
Plot 2D cell values
"""
function plot_cell_values(mesh::Mesh2D, cell_values::Vector{Float64}; title::String = "Finite Volume Plot")
    fig = Figure(fontsize = 18, size = (700, 600))
    ax  = Axis(fig[1, 1], title = title, aspect = DataAspect(), xlabel = "x", ylabel = "y")
    polys = [Polygon([Point2f(mesh.points[v]) for v in c]) for c in mesh.cells]
    m = poly!(ax, polys, color = cell_values, colormap = :viridis, strokewidth = 0.3, strokecolor = :black)
    Colorbar(fig[1, 2], m, label = "Cell Values")
    return fig
end

"""
Interactive animation of a time dependent problem on a 2D unstructured mesh.
Press 'm' to toggle mesh overlay.
"""
function animate_cell_values(
    mesh::Mesh2D,
    U_hist::Vector{Vector{Float64}};
    dt_hist = nothing
)
    t_hist = !isnothing(dt_hist) && length(dt_hist) == length(U_hist) ? cumsum(dt_hist) : nothing

    polys  = [Polygon([Point2f(mesh.points[v]) for v in c]) for c in mesh.cells]
    cmin   = minimum(minimum(u) for u in U_hist)
    cmax   = maximum(maximum(u) for u in U_hist)
    crange = cmin == cmax ? (cmin - 1.0, cmax + 1.0) : (cmin, cmax)

    fig = Figure(size = (750, 650))
    ax  = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
    sl  = Slider(fig[2, 1], range = 1:length(U_hist), startvalue = 1)

    obs = @lift(U_hist[$(sl.value)])
    m   = poly!(ax, polys, color = obs, colormap = :viridis, colorrange = crange, strokewidth = 0)
    Colorbar(fig[1, 2], m)

    mesh_visible = Observable(false)
    poly!(ax, polys, color = (:white, 0), strokecolor = :black, strokewidth = 0.5, visible = mesh_visible)

    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.m
            mesh_visible[] = !mesh_visible[]
        end
    end

    on(sl.value) do idx
        if !isnothing(t_hist)
            ax.title = "step $idx  |  dt = $(round(dt_hist[idx], sigdigits=3))  |  t = $(round(t_hist[idx], digits=4))"
        else
            ax.title = "Time Step: $(idx - 1)"
        end
    end

    GLMakie.activate!()
    screen = display(fig)
    @async for i in 1:length(U_hist)
        sl.value[] = i
        sleep(1/60)
    end
    wait(screen)
    return fig
end

"""
Save a 2D animation to a video file using CairoMakie
"""
function save_animation(
    mesh::Mesh2D,
    U_hist::Vector{Vector{Float64}},
    filename::String;
    dt_hist = nothing,
    framerate::Int = 60
)
    t_hist = !isnothing(dt_hist) && length(dt_hist) == length(U_hist) ? cumsum(dt_hist) : nothing

    polys  = [Polygon([Point2f(mesh.points[v]) for v in c]) for c in mesh.cells]
    cmin   = minimum(minimum(u) for u in U_hist)
    cmax   = maximum(maximum(u) for u in U_hist)
    crange = cmin == cmax ? (cmin - 1.0, cmax + 1.0) : (cmin, cmax)

    CairoMakie.activate!()

    fig = Figure(size = (700, 600))
    ax  = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x", ylabel = "y")
    obs = Observable(U_hist[1])
    m   = poly!(ax, polys, color = obs, colormap = :viridis, colorrange = crange, strokewidth = 0.3, strokecolor = :black)
    Colorbar(fig[1, 2], m)

    record(fig, filename, 1:length(U_hist); framerate) do idx
        obs[] = U_hist[idx]
        ax.title = !isnothing(t_hist) ?
            "step $idx  |  dt = $(round(dt_hist[idx], sigdigits=3))  |  t = $(round(t_hist[idx], digits=4))" :
            "Time Step: $(idx - 1)"
    end

    GLMakie.activate!()
    println("Saved animation to $filename")
end
