# solve the wave equation uₜₜ - c²uₓₓ = 0
# let v = uₜ and w = uₓ then we get a conservative system Uₜ + A*∇U = 0, with U = (v,w) and A = (0, -c^2; -1, 0)

# but it can be decomposed into another system : c²=κ/ρ  and U = (p,u) and A = (0, 1/ρ; κ, 0)

include("mesh1D.jl")
include("1D-animation.jl")

"""
computes the flux for the transport
"""
function get_flux_wave_1D(UL::Vector{Float64}, UR::Vector{Float64}, κ::Float64, ρ::Float64)
    p1, u1 = UL 
    p2, u2 = UR
    F1 = (0.5/ρ)*(u2+u1) - (0.5*κ/ρ)*(p2-p1)
    F2 = (0.5*κ)*(p2+p1) - (0.5*κ/ρ)*(u2-u1)
    return (F1, F2)
end

"""
Solves a wave step
"""
function step_wave!(new_cell_values, cell_values, mesh, dt, dx, step, κ::Float64, ρ::Float64; Ubr, Ubl)
    # interior nodes
    for (CL,CR) in mesh.Faces
        if CL!=0 && CR!=0
            UL = cell_values[CL]
            UR = cell_values[CR]
            F = get_flux_wave_1D(UL, UR, κ, ρ)
            new_cell_values[CL] -= (dt/dx) * F
            new_cell_values[CR] += (dt/dx) * F
        end
    end

    # boundary conditions
    # TODO: non periodic bcs
    if length(mesh.BoundaryFaces) > 0
       skip 
    end
end

"""
Solve the wave equation
"""
function solve_wave_1D(mesh::Mesh1D, n_timesteps::Int, dx::Real, dt::Real, u0::Function, p0::Function, κ::Real, ρ::Real)
    c = sqrt(κ/ρ)
    v0 = 0.5*p0-(0.5/(ρ*c))*u0
    w0 = 0.5*p0+(0.5/(ρ*c))*u0

    x0, x1 = mesh.x[1], mesh.x[end]
    xmid = (mesh.x[1:end-1] + mesh.x[2:end]) / 2
    
    u = u0.(xmid)
    u_new = copy(u)
    
    u_hist = [copy(u)]
    u_exact_hist = [copy(u)]

    p = p0.(xmid)
    p_new = copy(p)
    p_history = [copy(cell_values)]
    p_exact_history = [copy(cell_values)]

    for step in 1:n_timesteps
        step_wave!(new_cell_values, cell_values, mesh, dt, dx, step, κ, ρ)
        cell_values .= new_cell_values

        # exact solution 
        if length(mesh.BoundaryFaces) == 0
            vtrue = v0.(mod.(xmid .- c*dt*step, x1-x0))
            wtrue = w0.(mod.(xmid .+ c*dt*step, x1-x0))
            ptrue = vtrue.+wtrue
            utrue = (ρ*c)*(wtrue.-vtrue)
        end

        # save states
        push!(U_history, copy(cell_values))
        push!(U_exact_history, copy(utrue))
    end
    
    return xmid, U_history, U_exact_history

    

end




function main()
    # Parameters
    x0 = 0.0
    x1 = 1.0
    N = 100

    κ = 8.
    ρ = 2.
    c = sqrt(κ/ρ)

    # Initial condition
    p0(x) = 0.2 < x < 0.4 ? 100*exp(0.1/((x-0.2)*(x-0.4))) : 0
    u0(x) = -ρ*c*p0(x)

    # Generate mesh
    periodic = true
    mesh = generate_1DMesh(x0, x1, N, periodic)
    xmid = (mesh.x[1:end-1] + mesh.x[2:end]) / 2


    # display initial conditions
    display(plot1D(xmid, p0.(xmid), u0.(xmid), title="Initial condition"))
    # characteristic variables
    v0(x) = 0.5*p0(x) - u0(x)*0.5/(ρ*c)
    w0(x) = 0.5*p0(x) + u0(x)*0.5/(ρ*c)
    display(plot1D(xmid, v0.(xmid), w0.(xmid), title="Initial condition for the characteristic variables"))

    
    # solve the wave equation
    xmid, U_hist, U_exact_hist = solve_wave_1D(mesh, n_timesteps, dx, dt, u0, p0, κ, ρ)

end




main()
