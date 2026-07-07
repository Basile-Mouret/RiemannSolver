using FiniteVolumes
using StaticArrays

x0, x1 = 0.0, 1.5
N = 100

mesh = generate_1DMesh(x0, x1, N, false)
eq = Burgers1D()
bcs = Dict("left" => Outflow(), "right" => Outflow())

function u0(x)
    if x <= 0.5
        return SVector(-0.5)
    elseif x >= 1
        return SVector(0.0)
    else
        return SVector(1.0)
    end
end

CFL = 0.9
final_time = 0.5
max_time_steps = 1000

solve(mesh, eq, bcs, u0; max_time_steps = max_time_steps, CFL = CFL, final_time = final_time, output_dir = "out/burgers_Toro1D")
