using CairoMakie

"""
Plot 1D cell values, add true solution if available
"""
function plot_cell_values(mesh::Mesh1D, u::Vector{Float64}; u_exact = nothing, title::String = "Finite Volume Plot")
    fig = Figure(size = (800, 500))
    ax = Axis(fig[1, 1], title = title, xlabel = "x", ylabel = "U")
    stairs!(ax, mesh.cell_centers, u, step = :center, label = "Numerical", linewidth = 2)
    if !isnothing(u_exact)
        stairs!(ax, mesh.cell_centers, u_exact, step = :center, label = "Exact", linestyle = :dash, linewidth = 2)
    end
    axislegend(ax)
    return fig
end

"""
animates a time dependent problem on a 1D Mesh
"""
function animate_cell_values(
    mesh::Mesh1D,
    U_hist::Vector{Vector{Float64}},
    filename::String;
    U_exact_hist = nothing,
    dt_hist = nothing
)
    if !isnothing(dt_hist) && length(dt_hist) == length(U_hist)
        t_hist = cumsum(dt_hist)
    else
        t_hist = nothing
    end

    fig_anim = Figure(size = (800, 500))
    ax_anim = Axis(fig_anim[1, 1])
    obs_numerical = Observable(U_hist[1])
    stairs!(ax_anim, mesh.cell_centers, obs_numerical, step = :center, label = "Numerical")
    if !isnothing(U_exact_hist)
        obs_exact = Observable(U_exact_hist[1])
        stairs!(ax_anim,  mesh.cell_centers, obs_exact, step = :center, label = "Exact", linestyle = :dash)
    end
    axislegend(ax_anim)
    record(fig_anim, filename, 1:length(U_hist), framerate = 60) do frame_idx
        obs_numerical[] = U_hist[frame_idx]
        if !isnothing(U_exact_hist)
            obs_exact[] = U_exact_hist[frame_idx]
        end
        if !isnothing(dt_hist) && !isnothing(t_hist)
            t_val = t_hist[frame_idx]
            dt_val = dt_hist[frame_idx]
            title_str = "step $(frame_idx)  |  dt = $(round(dt_val, sigdigits=3))  |  t = $(round(t_val, digits=4))"
        else
            title_str = "Time Step: $(frame_idx - 1)"
        end
        ax_anim.title = title_str
    end
    println("Saved animation to $filename")
end
