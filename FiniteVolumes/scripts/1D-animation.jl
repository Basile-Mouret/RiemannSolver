using CairoMakie

"""
plot a 1D solution
"""
function plot1D(xmid::Vector{Float64}, u::Vector{Float64}; u_exact=nothing, title::String="Finite Volume Plot")
    # Initialize Figure
    fig = Figure(size = (800, 500))
    
    # Create the axis within the figure
    ax = Axis(fig[1, 1], title = title, xlabel="x", ylabel="U")
    
    # Plot the numerical solution
    stairs!(ax, xmid, u, step=:center, label="Numerical", linewidth=2)
    
    # Plot the exact solution if it was provided
    if !isnothing(u_exact)
        stairs!(ax, xmid, u_exact, step=:center, label="Exact", linestyle=:dash, linewidth=2)
    end
    
    # Add the legend to the axis
    axislegend(ax)
    
    # Return the figure
    return fig
end

function plot1D(xmid::AbstractVector{<:Real}, u::AbstractVector{<:Real}, v::AbstractVector{<:Real}; u_exact=nothing, v_exact=nothing, title::String="Finite Volume Plot")
    # Initialize Figure
    fig = Figure(size = (800, 1000))
    
    # Create the axis within the figure
    ax1 = Axis(fig[1, 1], title = title, xlabel="x", ylabel="U")
    ax2 = Axis(fig[2, 1], title = title, xlabel="x", ylabel="V")
    
    # Plot the numerical solution
    stairs!(ax1, xmid, u, step=:center, label="Numerical", linewidth=2)
    stairs!(ax2, xmid, v, step=:center, label="Numerical", linewidth=2)
    
    # Plot the exact solution if it was provided
    if !isnothing(u_exact)
        stairs!(ax1, xmid, u_exact, step=:center, label="Exact", linestyle=:dash, linewidth=2)
    end
    if !isnothing(v_exact)
        stairs!(ax2, xmid, v_exact, step=:center, label="Exact", linestyle=:dash, linewidth=2)
    end
    
    # Add the legend to the axis
    axislegend(ax1)
    axislegend(ax2)
    
    # Return the figure
    return fig
end

"""
animation function for 1D finite volumes using CairoMakie
is uses the stairs method to plot the piecewise constant data.
"""
function animate_1D_solution(xmid, U_hist, U_exact_hist, filename::String)
    fig_anim = Figure(size = (800, 500))
    # ax_anim = Axis(fig_anim[1, 1], limits = (minimum(xmid), maximum(xmid), -1.5, 1.5))
    ax_anim = Axis(fig_anim[1, 1])

    # Observables bound to the first frame
    obs_numerical = Observable(U_hist[1])
    obs_exact = Observable(U_exact_hist[1])

    stairs!(ax_anim, xmid, obs_numerical, step=:center, label="Numerical")
    stairs!(ax_anim, xmid, obs_exact, step=:center, label="Exact", linestyle=:dash)
    axislegend(ax_anim)

    record(fig_anim, filename, 1:length(U_hist), framerate=60) do frame_idx
        ax_anim.title = "Time Step: $(frame_idx - 1)"
        obs_numerical[] = U_hist[frame_idx]
        obs_exact[] = U_exact_hist[frame_idx]
    end
    
    println("Saved animation to $filename")
end
