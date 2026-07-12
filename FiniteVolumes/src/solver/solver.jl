"""
interface to solve a hyperbolic problem using finite volumes
"""
function solve(
    mesh::M,
    eq::E,
    boundary_conditions::B,
    ic::IC;
    max_time_steps::Int,
    final_time::Float64,
    CFL::Float64 = 0.9,
    output_dir::String = "out/simulation",
    dt_out::Float64 = final_time / 100,
    dt_max::Float64 = Inf,
    Verbose::Bool=true,
    n_info::Int=100
) where {M<:AbstractMesh, E<:AbstractEquation, B, IC}

    nvars = num_vars(eq)

    N = length(mesh.cells)
    values = Matrix{Float64}(undef, N, nvars)
    for i in 1:N
        values[i, :] .= ic(mesh.cell_centers[i])
    end
    new_values = copy(values)

    step = 0
    t = 0.0

    outdir   = rstrip(output_dir, '/')
    filename = basename(outdir)
    writer = VTKStreamWriter(mesh, eq; outdir=outdir, filename=filename, dt_out=dt_out)
    write_frame!(writer, values, t)

    if Verbose
        println("Solve started")
        solve_start_time = time()
    end
    while t < final_time && step < max_time_steps
        dt = min(dt_max, compute_dt(mesh, eq, values, CFL))
        if t + dt > final_time
            dt = final_time - t
        end

        new_values .= values

        explicit_euler_step!(new_values, values, mesh, eq, boundary_conditions, dt, t)
        values .= new_values


        if Verbose && step % (max_time_steps ÷ n_info) == 1
            _print_sim_info(step, max_time_steps, t, final_time, solve_start_time)
        end

        t += dt
        step += 1

        maybe_write!(writer, values, t)
    end

    close_writer!(writer)
    if Verbose
        comp_time = time()-solve_start_time
        println("Solve finished in ",comp_time,"s in ", step, " steps.")
    end
end

"""
prints simulation information during computation
"""
function _print_sim_info(step::Int, max_time_steps::Int, t::Float64, final_time::Float64, solve_start_time::Float64)
    elapsed = time() - solve_start_time
    elapsed_str = convert_time(elapsed)
    progress = step / max_time_steps

    if step == 0 || t == 0
        eta_str = "N/A"
    else
        eta = elapsed * min((final_time - t) / t, (1 - progress) / progress)
        eta_str = convert_time(eta)
    end

    println("step ", step, "/", max_time_steps,
             " (", round(100*progress, digits=1), "%), ",
             "t = ", round(t, digits=4), "/", round(final_time, digits=4),
             " (", round(100*t/final_time, digits=1), "%), ",
             "elapsed = ", elapsed_str, "s | ",
             "eta = ", eta_str)
end

"""
converts seconds into a string in the hh:mm:ss format
"""
function convert_time(seconds::Float64)
    h = floor(Int,seconds / 3600)
    m = floor(Int,(seconds % 3600) / 60)
    s = floor(Int,seconds % 60)
    return string(h, ":", lpad(m, 2, '0'), ":", lpad(s, 2, '0'))
end

