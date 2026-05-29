using CairoMakie

"""
ploting function for 1D finite volumes using CairoMakie
is uses the stairs method to plot the piecewise constant data.
"""
function animate_1D_solution(xmid, U_hist, U_exact_hist, filename::String)
    fig_anim = Figure(size = (800, 500))
    ax_anim = Axis(fig_anim[1, 1], limits = (minimum(xmid), maximum(xmid), -1.5, 1.5))

    # Observables bound to the first frame
    obs_numerical = Observable(U_hist[1])
    obs_exact = Observable(U_exact_hist[1])

    stairs!(ax_anim, xmid, obs_numerical, step=:post, label="Numerical")
    stairs!(ax_anim, xmid, obs_exact, step=:post, label="Exact", linestyle=:dash)
    axislegend(ax_anim)

    record(fig_anim, filename, 1:length(U_hist), framerate=60) do frame_idx
        ax_anim.title = "Time Step: $(frame_idx - 1)"
        obs_numerical[] = U_hist[frame_idx]
        obs_exact[] = U_exact_hist[frame_idx]
    end
    
    println("Saved animation to $filename")
end
