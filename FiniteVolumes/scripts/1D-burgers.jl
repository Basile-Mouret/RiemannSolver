
include("mesh1D.jl")
include("1D-animation.jl")

"""
computes the flux for the transport
"""
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

"""
Solves a transport step
"""
function step_transport!(newCellValues, CellValues, mesh, c, dt, dx, step, flux_type, ub)
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
            F_in = c * ub(dt*step)
            F_out = c * CellValues[end]
            newCellValues[1] += (dt/dx) * F_in
            newCellValues[end] -= (dt/dx) * F_out
        else
            F_in = c * ub(dt*step)
            F_out = c * CellValues[1]
            newCellValues[end] -= (dt/dx) * F_in
            newCellValues[1] += (dt/dx) * F_out
        end
    end
end

"""
Solves transport over some timesteps
"""
function solve_transport_1D(mesh::Mesh1D, c::Float64, n_timesteps::Int, dx::Float64, dt::Float64, u0::Function; ub::Function, flux_type::String="upwind")
    x0, x1 = mesh.x[1], mesh.x[end]
    xmid = (mesh.x[1:end-1] + mesh.x[2:end]) / 2
    
    CellValues = u0.(xmid)
    newCellValues = copy(CellValues)
    
    println("Initial L2 Error: ", compute_L2_1D(mesh, CellValues, u0))

    U_history = [copy(CellValues)]
    U_exact_history = [copy(CellValues)]

    for step in 1:n_timesteps
        step_transport!(newCellValues, CellValues, mesh, c, dt, dx, step, flux_type, ub)
        CellValues .= newCellValues

        # exact solution 
        if length(mesh.BoundaryFaces) == 0
            utrue = u0.(mod.(xmid .- c*dt*step, x1-x0))
        else
            if c > 0
                utrue = [x <= c*dt*step ? ub(dt*step - x/c) : u0(x - c*dt*step) for x in xmid]
            else
                utrue = [x >= x1 + c*dt*step ? ub(dt*step + (x1-x)/c) : u0(x - c*dt*step) for x in xmid]
            end
        end

        # save states
        push!(U_history, copy(CellValues))
        push!(U_exact_history, copy(utrue))
    end
    
    return xmid, U_history, U_exact_history
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

    # solve
    xmid, U_hist, U_exact_hist = solve_transport_1D(mesh, c, n_timesteps, dx, dt, u0, ub=ub, flux_type="upwind")

    # plot initial condition
    display(plot1D(xmid, u0.(xmid); title="Initial condition"))

    # animate
    anim_file = "media/advection_1d.mp4"
    animate_1D_solution(xmid, U_hist, U_exact_hist, anim_file)
    
    # open animation
    run(`xdg-open $anim_file`)
end

main()
