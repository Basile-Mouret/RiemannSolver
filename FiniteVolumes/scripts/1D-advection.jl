using Plots
gr()


"""
1D Mesh
"""
struct Mesh1D
    x::Vector{Float64}
    Cells::Vector{Tuple{Int, Int}}
    Faces::Vector{Tuple{Int, Int}}
    BoundaryFaces::Vector{Int}
end

"""
generates a 1D Mesh
"""
function generate_1DMesh(x0::Float64, x1::Float64, N::Int, periodic::Bool)
    x = range(x0, x1, length=N+1)
    Cells = [(i,i+1) for i in 1:N]
    Faces = [(i,i+1) for i in 1:N-1]
    if periodic
        push!(Faces, (N, 1))
        BoundaryFaces = []
    else
        BoundaryFaces = [1, N]
    end
    return Mesh1D(x, Cells, Faces, BoundaryFaces)
end

"""
computes a 1D, second order quadrature
"""
function quadrature_1D(x1::Float64, x2::Float64, f::Function)
    return (x2 - x1)/2 * (f(-(x2-x1)/(2*sqrt(3)) + (x1+x2)/2) + f((x2-x1)/(2*sqrt(3)) + (x1+x2)/2))
end

"""
computes L2 error on a 1D mesh
"""
function compute_L2_1D(mesh::Mesh1D, u::Vector{Float64}, u_exact::Function)
    error = 0.0
    for (i, (n1, n2)) in enumerate(mesh.Cells)
        x1 = mesh.x[n1]
        x2 = mesh.x[n2]
        error += quadrature_1D(x1, x2, x -> (u[i] - u_exact(x))^2)
    end
    return sqrt(error)
end

function fluxTransport1D(c, UL, UR)
    return c>0 ? c*UL : c*UR
end

function solve_transport_1D(x0::Float64, x1::Float64, N::Int, periodic::Bool)
    # Generate mesh
    mesh = generate_1DMesh(x0, x1, N, periodic)
    # Initial condition
    u0(x) = sin(2π*x)
    c = 1.5
    dx = (x1-x0)/N
    CFL = 0.8
    dt = CFL*dx/abs(c)
    n_timesteps = 100
    T = dt*n_timesteps  # Final simulation time

    xmid = (mesh.x[1:end-1] + mesh.x[2:end]) / 2
    CellValues = u0.(xmid)

    plt_init = plot(xmid, u0.(xmid), xlims=(x0, x1), ylims=(-1.5, 1.5), label="Initial Condition")
    
    plot!(plt_init, mesh.x[1:end-1], CellValues, seriestype = :steppost, label="Cell Values")
    plot!(plt_init, mesh.x, u0.(mesh.x), label="Exact Solution", linestyle=:dash)
    
    display(plt_init)

    println("L2 Error: ", compute_L2_1D(mesh, CellValues, u0))

    newCellValues = copy(CellValues)
    anim = @animate for step in 1:n_timesteps
        #update the cells
        # this can be computed in parallel
        for (CL,CR) in mesh.Faces
            # solve Riemann Problem
            UL = CellValues[CL]
            UR = CellValues[CR]
            F = fluxTransport1D(c, UL, UR)
            newCellValues[CL] -= (dt/dx) * F
            newCellValues[CR] += (dt/dx) * F
        end

        CellValues .= newCellValues
        utrue = u0.(mod.(xmid.-c*step*dt,1.0))
        plot(xmid, CellValues,
             xlims=(x0, x1), 
             ylims=(-1.5, 1.5), # Adjusted to fit the sine wave
             seriestype = :steppost,
             label="Numerical")

        plot!(xmid, utrue, 
             xlims=(x0, x1), 
             ylims=(-1.5, 1.5),
             seriestype = :steppost,
             label="Exact")
    end
    
    # Actually generate the file!
    gif(anim, "advection_1d.gif", fps=60)
    println("Saved animation to advection_1d.gif")

end


function main()
    # Parameters
    x0 = 0.0
    x1 = 1.0
    N = 100
    periodic = true

    solve_transport_1D(x0,x1, N, periodic)
end

main()

