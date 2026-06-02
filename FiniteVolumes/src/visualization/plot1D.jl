using CairoMakie

"""
Plot 1D cell values, add true solution if available
"""
function plot1D(xmid::Vector{Float64}, u::Vector{Float64}; u_exact = nothing, title::String = "Finite Volume Plot")
    fig = Figure(size = (800, 500))
    ax = Axis(fig[1, 1], title = title, xlabel = "x", ylabel = "U")
    stairs!(ax, xmid, u, step = :center, label = "Numerical", linewidth = 2)
    if !isnothing(u_exact)
        stairs!(ax, xmid, u_exact, step = :center, label = "Exact", linestyle = :dash, linewidth = 2)
    end
    axislegend(ax)
    return fig
end

"""
Plot 1D cell values for a system, add true solution if available
"""
function plot1D(
    xmid::AbstractVector{<:Real},
    u::AbstractVector{<:Real},
    v::AbstractVector{<:Real};
    u_exact = nothing,
    v_exact = nothing,
    title::String = "Finite Volume Plot"
)
    fig = Figure(size = (800, 1000))
    ax1 = Axis(fig[1, 1], title = title, xlabel = "x", ylabel = "U")
    ax2 = Axis(fig[2, 1], title = title, xlabel = "x", ylabel = "V")
    stairs!(ax1, xmid, u, step = :center, label = "Numerical", linewidth = 2)
    stairs!(ax2, xmid, v, step = :center, label = "Numerical", linewidth = 2)
    if !isnothing(u_exact)
        stairs!(ax1, xmid, u_exact, step = :center, label = "Exact", linestyle = :dash, linewidth = 2)
    end
    if !isnothing(v_exact)
        stairs!(ax2, xmid, v_exact, step = :center, label = "Exact", linestyle = :dash, linewidth = 2)
    end
    axislegend(ax1)
    axislegend(ax2)
    return fig
end

"""
animates a time dependent problem on a 1D Mesh
"""
function animate_1D_solution(
    xmid::Vector{Float64},
    U_hist::Vector{Vector{Float64}},
    filename::String;
    U_exact_hist = nothing
)
    fig_anim = Figure(size = (800, 500))
    ax_anim = Axis(fig_anim[1, 1])
    obs_numerical = Observable(U_hist[1])
    stairs!(ax_anim, xmid, obs_numerical, step = :center, label = "Numerical")
    if !isnothing(U_exact_hist)
        obs_exact = Observable(U_exact_hist[1])
        stairs!(ax_anim, xmid, obs_exact, step = :center, label = "Exact", linestyle = :dash)
    end
    axislegend(ax_anim)
    record(fig_anim, filename, 1:length(U_hist), framerate = 60) do frame_idx
        ax_anim.title = "Time Step: $(frame_idx - 1)"
        obs_numerical[] = U_hist[frame_idx]
        if !isnothing(U_exact_hist)
            obs_exact[] = U_exact_hist[frame_idx]
        end
    end
    println("Saved animation to $filename")
end
