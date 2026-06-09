"""
Interface to solve a hyperbolic problem using finite volumes
"""
function solve(
    mesh::Mesh1D,
    eq::AbstractEquation1D,
    bcs,
    ic::Function;
    max_time_steps::Int,
    final_time::Float64,
    CFL::Float64 = 0.9,
    compute_exact::Bool = false
)
    x0, x1 = mesh.points[1], mesh.points[end]
    nvars = num_vars(eq)

    N = length(mesh.points) - 1
    values = Matrix{Float64}(undef, N, nvars)
    for i in 1:N
        values[i, :] .= ic(mesh.cells_center[i])
    end
    new_values = copy(values)

    U_history = [copy(values)]
    U_exact_hist = [copy(values)]
    dt_hist = Float64[]

    step = 0
    t = 0.0
    while t < final_time && step < max_time_steps
        dt = CFL * mesh.cells_length[1] / max_wave_speed(eq, values, mesh) # to generalize with a helper function for unregular grids (compute cfl)
        if t + dt > final_time
            dt = final_time - t
        end

        new_values .= values

        explicit_euler_step!(new_values, values, mesh, eq, bcs, dt, t)
        values .= new_values
        t += dt
        step += 1

        push!(dt_hist, dt)
        push!(U_history, copy(values))

        if compute_exact
            exact = zeros(N, nvars)
            exact_solution!(exact, eq, mesh.cells_center, ic, bcs, x0, x1, t)
            push!(U_exact_hist, copy(exact))
        else
            push!(U_exact_hist, copy(values))
        end
    end

    return U_history, U_exact_hist, [0.0; dt_hist]
end

