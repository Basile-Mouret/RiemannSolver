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
    compute_exact::Bool = true
)
    xmid = cell_centers(mesh)
    dx = cell_width(mesh)
    x0, x1 = mesh.x[1], mesh.x[end]
    nvars = num_vars(eq)

    N = length(mesh.x) - 1
    values = Matrix{Float64}(undef, N, nvars)
    for i in 1:N
        values[i, :] .= ic(xmid[i])
    end
    new_values = copy(values)

    U_history = [copy(values)]
    U_exact_hist = [copy(values)]

    step = 0
    t = 0.0
    while t < final_time && step < max_time_steps
        # compute adaptive dt
        dt = CFL*dx/max_wave_speed(eq, values, mesh)
        t = dt * (step - 1)

        new_values .= values

        explicit_euler_step!(new_values, values, mesh, eq, bcs, dt, dx, t)
        values .= new_values

        push!(U_history, copy(values))

        if compute_exact
            t_end = dt * step
            exact = zeros(N, nvars)
            exact_solution!(exact, eq, xmid, ic, bcs, x0, x1, t_end)
            push!(U_exact_hist, copy(exact))
        else
            push!(U_exact_hist, copy(values))
        end
        t+=dt
        step+=1
    end

    return xmid, U_history, U_exact_hist
end
