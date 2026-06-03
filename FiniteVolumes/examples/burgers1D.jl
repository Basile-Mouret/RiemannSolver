using CairoMakie
using FiniteVolumes

x0, x1 = 0.0, 1.0
N = 1000

mesh = generate_1DMesh(x0, x1, N, false)
eq = Burgers1D(1.)
bcs = Dict(:left => Outflow(), :right => Outflow())

function u0(x)
    return sin(2*π*x)
end

CFL = 0.4

xmid, U_hist, U_exact_hist = solve(mesh, eq, bcs, u0; nsteps = 100, CFL = CFL)

u0_vals = [u0(x)[1] for x in xmid]
display(plot1D(xmid, u0_vals; title = "Initial condition"))

anim_file = "media/burgers_1d.mp4"
U_scalar = [mat[:, 1] for mat in U_hist]
U_exact_scalar = [mat[:, 1] for mat in U_exact_hist]
animate_1D_solution(xmid, U_scalar, anim_file; U_exact_hist = U_exact_scalar)

run(`xdg-open $(anim_file)`)
