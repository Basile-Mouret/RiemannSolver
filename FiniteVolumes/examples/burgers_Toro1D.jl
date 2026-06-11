using CairoMakie
using FiniteVolumes

x0, x1 = 0.0, 1.5
N = 100

mesh = generate_1DMesh(x0, x1, N, false)
eq = Burgers1D()
bcs = Dict("left" => Outflow(), "right" => Outflow())

function u0(x)
    if x<= 0.5
        return -0.5
    elseif x>=1
        return 0.0
    else
        return 1.0
    end
end

CFL = 0.9
final_time = 0.5
max_time_steps = 1000

U_hist, dt_hist = solve(mesh, eq, bcs, u0; max_time_steps = max_time_steps, CFL = CFL, final_time=final_time)

u0_vals = [u0(x)[1] for x in mesh.cell_centers]
display(plot_cell_values(mesh, u0_vals; title = "Initial condition"))

anim_file = "media/burgers_1d.mp4"
U_scalar = [mat[:, 1] for mat in U_hist]
animate_cell_values(mesh, U_scalar, anim_file; dt_hist = dt_hist)

run(`xdg-open $(anim_file)`, wait=false)
