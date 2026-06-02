# solve the wave equation uₜₜ - c²uₓₓ = 0
# let v = uₜ and w = uₓ then we get a conservative system Uₜ + A*∇U = 0, with U = (v,w) and A = (0, -c^2; -1, 0)

# but it can be decomposed into another system : c²=κ/ρ  and U = (p,u) and A = (0, 1/ρ; κ, 0)

include("mesh1D.jl")
include("1D-animation.jl")

"""
computes the flux for the transport
"""
function compute_flux_wave_1D(pl::T, ul::T, pr::T, ur::T, κ::T, ρ::T) where {T<:Real}
    F_p = (0.5 / ρ) * (ur + ul) - (0.5 * sqrt(κ/ρ)) * (pr - pl)
    F_u = (0.5 * κ) * (pr + pl) - (0.5 * sqrt(κ/ρ)) * (ur - ul)
    return F_p, F_u
end

"""
Solves a wave step
"""
function step_wave!(new_cell_values::Matrix{T}, cell_values::Matrix{T}, mesh::Mesh1D, dt::T, dx::T, step::Int, κ::T, ρ::T; Ubr::Function=(x->0), Ubl::Function=(x->0)) where {T<:Real}
    # interior nodes
    for (CL,CR) in mesh.Faces
        if CL!=0 && CR!=0
            pl = cell_values[CL, 1]
            ul = cell_values[CL, 2]

            pr = cell_values[CR, 1]
            ur = cell_values[CR, 2]
        elseif CL == 0
            pl = cell_values[CR, 1]
            ul = -cell_values[CR, 2]

            pr = cell_values[CR, 1]
            ur = cell_values[CR, 2]
        elseif CR == 0
            pl = cell_values[CL, 1]
            ul = cell_values[CL, 2]

            pr = cell_values[CL, 1]
            ur = -cell_values[CL, 2]
        end

        F_p, F_u = compute_flux_wave_1D(pl, ul, pr, ur, κ, ρ)

        if CL != 0
            new_cell_values[CL, 1] -= (dt/dx) * F_p
            new_cell_values[CL, 2] -= (dt/dx) * F_u
        end

        if CR != 0
            new_cell_values[CR, 1] += (dt/dx) * F_p
            new_cell_values[CR, 2] += (dt/dx) * F_u
        end
    end
end

"""
Solve the wave equation
"""
function solve_wave_1D(mesh::Mesh1D, n_timesteps::Int, dx::T, dt::T, u0::Function, p0::Function, κ::T, ρ::T) where {T<:Real}
    c = sqrt(κ/ρ)
    v0(x) = 0.5*p0(x)-(0.5/(ρ*c))*u0(x)
    w0(x) = 0.5*p0(x)+(0.5/(ρ*c))*u0(x)

    x0, x1 = mesh.x[1], mesh.x[end]
    xmid = (mesh.x[1:end-1] + mesh.x[2:end]) / 2
    
    CellValues = Float64.(hcat(p0.(xmid), u0.(xmid)))
    newCellValues = copy(CellValues)
    
    U_history = [copy(CellValues)]
    U_exact_history = [copy(CellValues)]

    for step in 1:n_timesteps
        step_wave!(newCellValues, CellValues, mesh, dt, dx, step, κ, ρ)
        CellValues .= newCellValues

        # exact solution 
        if length(mesh.BoundaryFaces) == 0
            vtrue = v0.(mod.(xmid .+ c*dt*step, x1-x0))
            wtrue = w0.(mod.(xmid .- c*dt*step, x1-x0))
        else
            vtrue = 0 .*xmid
            wtrue = 0 .*xmid
        end
        U_true = hcat(vtrue.+wtrue, (ρ*c)*(wtrue.-vtrue))

        # save states
        push!(U_history, copy(CellValues))
        push!(U_exact_history, copy(U_true))
    end
    
    return xmid, U_history, U_exact_history
end

function compute_entropy(cell_values::Matrix{T}, dx::T, ρ::T, κ::T) where {T<:Real}
    p = cell_values[:, 1]
    u = cell_values[:, 2]
    return 0.5*dx*(sum(p.*p+ (u.*u)/(ρ*κ))) 
end

function main()
    # Parameters
    x0 = 0.0
    x1 = 1.0
    N = 100

    κ = 8.
    ρ = 2.
    c = sqrt(κ/ρ)

    dx = (x1-x0)/N
    CFL = 0.8
    dt = CFL*dx/abs(c)
    n_timesteps = 1000


    # Initial condition
    p0(x) = 0.2 < x < 0.4 ? 100*exp(0.1/((x-0.2)*(x-0.4))) : 0
    u0(x) = 0*x  #-ρ*c*p0(x)

    # Generate mesh
    periodic = false
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

    # entropy graph over time
    entropy_hist = Float64[]
    for mat in U_hist
        entropy = compute_entropy(mat, dx, ρ, κ)
        push!(entropy_hist, entropy)
    end
    
    display(plot(entropy_hist))



    # animation
    
    p_hist = [mat[:, 1] for mat in U_hist]
    u_hist = [mat[:, 2] for mat in U_hist]

    p_exact_hist = [mat[:, 1] for mat in U_exact_hist]
    u_exact_hist = [mat[:, 2] for mat in U_exact_hist]

    anim_file = "media/wave_u_1d.mp4"
    animate_1D_solution(xmid, u_hist, u_exact_hist, anim_file)
    run(`xdg-open $anim_file`)

    anim_file = "media/wave_p_1d.mp4"
    animate_1D_solution(xmid, p_hist, p_exact_hist, anim_file)
    run(`xdg-open $anim_file`)



end


main()

