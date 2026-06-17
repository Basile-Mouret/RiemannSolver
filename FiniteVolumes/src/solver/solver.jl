"""
Interface to solve a hyperbolic problem using finite volumes
"""
function solve(
    mesh::M,
    eq::E,
    bcs::B,
    ic::IC;
    max_time_steps::Int,
    final_time::Float64,
    CFL::Float64 = 0.9,
    output_dir::String = "out/out",
    dt_out::Float64 = final_time / 100,
) where {M<:AbstractMesh, E<:AbstractEquation, B, IC}

    nvars = num_vars(eq)

    N = length(mesh.cells)
    values = Matrix{Float64}(undef, N, nvars)
    for i in 1:N
        values[i, :] .= ic(mesh.cell_centers[i])
    end
    new_values = copy(values)

    #U_history = [copy(values)]
    #timesteps = Float64[]

    step = 0
    t = 0.0

    outdir   = dirname(output_dir)
    filename = basename(output_dir)
    isempty(outdir) && (outdir = ".")
    writer = VTKStreamWriter(mesh, eq; outdir=outdir, filename=filename, dt_out=dt_out)
    write_frame!(writer, values, t)

    while t < final_time && step < max_time_steps
        dt = compute_dt(mesh, eq, values, CFL)
        if t + dt > final_time
            dt = final_time - t
        end

        new_values .= values

        explicit_euler_step!(new_values, values, mesh, eq, bcs, dt, t)
        values .= new_values
        t += dt
        step += 1

        maybe_write!(writer, values, t)

        # push!(dt_hist, dt)
        # push!(U_history, copy(values))
    end

    close_writer!(writer)
end

