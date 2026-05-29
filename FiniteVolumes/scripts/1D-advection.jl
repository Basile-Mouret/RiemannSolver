using CairoMakie
include("mesh1D.jl")

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

function fluxTransport1D(c, UL, UR ; type="upwind")
    if type == "upwind"
        return c > 0 ? c*UL : c*UR
    elseif type == "centered"
        return 0.5*c*(UL+UR) 
    elseif type == "left"
        return c*UL 
    elseif type == "right"
        return c*UR 
    end
end

function solve_transport_1D(mesh::Mesh1D, c::Float64, n_timesteps::Int, dx::Float64, dt::Float64, u0::Function; ub::Function,  flux_type::String="upwind", animation_title::String="advection_1d.mp4")
    x0, x1 = mesh.x[1], mesh.x[end]

    xmid = (mesh.x[1:end-1] + mesh.x[2:end]) / 2
    CellValues = u0.(xmid)

    # initial plot
    fig_init = Figure(size = (800, 500))
    ax_init = Axis(fig_init[1, 1], limits = (x0, x1, -1.5, 1.5), title = "Initial State")
    
    lines!(ax_init, xmid, u0.(xmid), label="Initial Condition")
    stairs!(ax_init, mesh.x[1:end-1], CellValues, step=:post, label="Cell Values")
    lines!(ax_init, mesh.x, u0.(mesh.x), label="Exact Solution", linestyle=:dash)
    
    axislegend(ax_init)
    display(fig_init)

    println("L2 Error: ", compute_L2_1D(mesh, CellValues, u0))

    # animation
    newCellValues = copy(CellValues)
    
    # Create the animation figure
    fig_anim = Figure(size = (800, 500))
    ax_anim = Axis(fig_anim[1, 1], limits = (x0, x1, -1.5, 1.5), title = "Time Step: 0")

    # Wrap the changing arrays in Observables
    obs_numerical = Observable(CellValues)
    obs_exact = Observable(u0.(xmid))

    # Plot the observables once
    stairs!(ax_anim, xmid, obs_numerical, step=:post, label="Numerical")
    stairs!(ax_anim, xmid, obs_exact, step=:post, label="Exact", linestyle=:dash)
    axislegend(ax_anim)

    record(fig_anim, animation_title, 1:n_timesteps, framerate=60) do step
        
        # interior nodes
        for (CL,CR) in mesh.Faces
            if CL!=0 && CR!=0
                UL = CellValues[CL]
                UR = CellValues[CR]
                F = fluxTransport1D(c, UL, UR, type=flux_type)
                newCellValues[CL] -= (dt/dx) * F
                newCellValues[CR] += (dt/dx) * F
            end
        end

        # boundary conditions
        if length(mesh.BoundaryFaces) > 0
            if c > 0
                # flow from left to right.
                F_in = c * ub(dt*step)
                F_out = c * CellValues[end]
                
                newCellValues[1] += (dt/dx) * F_in
                newCellValues[end] -= (dt/dx) * F_out
            else
                # flow from right to left.
                F_in = c * ub(dt*step)
                F_out = c * CellValues[1]
                
                newCellValues[end] -= (dt/dx) * F_in
                newCellValues[1] += (dt/dx) * F_out
            end
        end
        CellValues .= newCellValues

        # compute the true analytic solution u(x,t) = u0(x-ct)
        if length(mesh.BoundaryFaces) == 0
            # periodic
            utrue = u0.(mod.(xmid .- c*dt*step, x1-x0))
        else
            # non-periodic
            if c > 0
                utrue = [x <= c*dt*step ? ub(dt*step - x/c) : u0(x - c*dt*step) for x in xmid]
            else
                utrue = [x >= x1 + c*dt*step ? ub(dt*step + (x1-x)/c) : u0(x - c*dt*step) for x in xmid]
            end
        end

        ax_anim.title = "Time Step: $step"

        obs_numerical[] = CellValues
        obs_exact[] = utrue
    end
    
    println("Saved animation to $animation_title")
end

function main()
    # Parameters
    x0 = 0.0
    x1 = 1.0
    N = 100

    # Initial condition
    #u0(x) = sin(2π*x)
    u0(x) = 0*x 
    # Boundary Conditions
    ub(t) = sin(4*π*t) 

    # propagation speed
    c = 1.5

    # Generate mesh
    periodic = false
    mesh = generate_1DMesh(x0, x1, N, periodic)

    dx = (x1-x0)/N
    CFL = 0.8
    dt = CFL*dx/abs(c)
    n_timesteps = 100


    solve_transport_1D(mesh, c, n_timesteps, dx, dt, u0, ub = ub,  flux_type="upwind", animation_title="advection_1d.mp4")
    run(`xdg-open advection_1d.mp4`)
end

main()
