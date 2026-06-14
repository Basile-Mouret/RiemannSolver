using GLMakie
import CairoMakie
GLMakie.activate!()

"""
Plot 1D cell values, add true solution if available
"""
function plot_cell_values(mesh::Mesh1D, u::Vector{Float64}; u_exact = nothing, title::String = "Finite Volume Plot")
    fig = Figure(size = (800, 500))
    ax = Axis(fig[1, 1], title = title, xlabel = "x", ylabel = "U")
    stairs!(ax, mesh.cell_centers, u, step = :center, label = "Numerical", linewidth = 2)
    all_vals = u
    if !isnothing(u_exact)
        stairs!(ax, mesh.cell_centers, u_exact, step = :center, label = "Exact", linestyle = :dash, linewidth = 2)
        all_vals = vcat(all_vals, u_exact)
    end
    axislegend(ax)
    ymin, ymax = extrema(all_vals)
    pad = 0.05 * max(ymax - ymin, eps())
    ylims!(ax, ymin - pad, ymax + pad)
    return fig
end

"""
Interactive animation of a time dependent problem on a 1D mesh
"""
function animate_cell_values(
    mesh::Mesh1D,
    U_hist::Vector{Vector{Float64}};
    U_exact_hist = nothing,
    dt_hist = nothing
)
    t_hist = !isnothing(dt_hist) && length(dt_hist) == length(U_hist) ? cumsum(dt_hist) : nothing

    all_vals = reduce(vcat, U_hist)
    !isnothing(U_exact_hist) && append!(all_vals, reduce(vcat, U_exact_hist))
    ymin, ymax = extrema(all_vals)
    pad = 0.05 * max(ymax - ymin, eps())

    fig = Figure(size = (700, 750))
    ax  = Axis(fig[1, 1], xlabel = "x", ylabel = "U", aspect = AxisAspect(1),
               limits = (nothing, nothing, ymin - pad, ymax + pad))
    sl  = Slider(fig[2, 1], range = 1:length(U_hist), startvalue = 1)

    obs_numerical = @lift(U_hist[$(sl.value)])
    stairs!(ax, mesh.cell_centers, obs_numerical, step = :center, label = "Numerical")

    if !isnothing(U_exact_hist)
        obs_exact = @lift(U_exact_hist[$(sl.value)])
        stairs!(ax, mesh.cell_centers, obs_exact, step = :center, label = "Exact", linestyle = :dash)
    end
    axislegend(ax)

    on(sl.value) do idx
        if !isnothing(t_hist)
            ax.title = "step $idx  |  dt = $(round(dt_hist[idx], sigdigits=3))  |  t = $(round(t_hist[idx], digits=4))"
        else
            ax.title = "Time Step: $(idx - 1)"
        end
    end

    mesh_visible = Observable(false)
    vlines!(ax, mesh.face_centers, color = :black, linewidth = 0.5, visible = mesh_visible)

    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press && event.key == Keyboard.m
            mesh_visible[] = !mesh_visible[]
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
Save a 1D animation to a video file using CairoMakie
"""
function save_animation(
    mesh::Mesh1D,
    U_hist::Vector{Vector{Float64}},
    filename::String;
    U_exact_hist = nothing,
    dt_hist = nothing,
    framerate::Int = 60
)
    t_hist = !isnothing(dt_hist) && length(dt_hist) == length(U_hist) ? cumsum(dt_hist) : nothing

    all_vals = reduce(vcat, U_hist)
    !isnothing(U_exact_hist) && append!(all_vals, reduce(vcat, U_exact_hist))
    ymin, ymax = extrema(all_vals)
    pad = 0.05 * max(ymax - ymin, eps())

    CairoMakie.activate!()

    fig = Figure(size = (800, 500))
    ax  = Axis(fig[1, 1], xlabel = "x", ylabel = "U",
               limits = (nothing, nothing, ymin - pad, ymax + pad))
    obs = Observable(U_hist[1])
    stairs!(ax, mesh.cell_centers, obs, step = :center, label = "Numerical")
    if !isnothing(U_exact_hist)
        obs_exact = Observable(U_exact_hist[1])
        stairs!(ax, mesh.cell_centers, obs_exact, step = :center, label = "Exact", linestyle = :dash)
    end
    axislegend(ax)

    record(fig, filename, 1:length(U_hist); framerate) do idx
        obs[] = U_hist[idx]
        !isnothing(U_exact_hist) && (obs_exact[] = U_exact_hist[idx])
        ax.title = !isnothing(t_hist) ?
            "step $idx  |  dt = $(round(dt_hist[idx], sigdigits=3))  |  t = $(round(t_hist[idx], digits=4))" :
            "Time Step: $(idx - 1)"
    end

    GLMakie.activate!()
    println("Saved animation to $filename")
end
