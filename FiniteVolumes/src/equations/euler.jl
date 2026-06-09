struct Euler1D <: AbstractEquation1D
    gamma::Float64
    solver::Symbol
end

function _cons_to_prim_euler_1D_ideal_gas(U::Vector{T}, gamma::T) where {T<:Real}
    rho = U[1]
    rhou = U[2]
    E = U[3]
    u = rhou / rho
    p = (gamma - 1.0) * (E - 0.5 * rhou^2 / rho)
    return [rho, u, p]
end

function _prim_to_cons_euler_1D_ideal_gas(W::Vector{T}, gamma::T) where {T<:Real}
    rho = W[1]
    u = W[2]
    p = W[3]
    rhou = rho * u
    E = p / (gamma - 1.0) + 0.5 * rho * u^2
    return [rho, rhou, E]
end

num_vars(::Euler1D) = 3

function max_wave_speed(eq::Euler1D, U::Matrix{Float64}, Mesh::Mesh1D )
s_max = 0.0
    
    for i in axes(U, 1)
        U_cell = U[i, :]
        
        rho, u, p = _cons_to_prim_euler_1D_ideal_gas(U_cell, eq.gamma)
        
        a = sqrt(eq.gamma * p / rho)
        s_max = max(s_max, abs(u) + a)
    end
    
    return s_max
end

function flux(eq::Euler1D, UL::Vector{Float64}, UR::Vector{Float64})
    
    WL = _cons_to_prim_euler_1D_ideal_gas(UL, eq.gamma)
    WR = _cons_to_prim_euler_1D_ideal_gas(UR, eq.gamma)
    
    # solve the Riemann problem at the interface (x/t = 0.0)
    if eq.solver == :exact
        rho, u, p = solve_riemann_exact(0.0, WL, WR, eq.gamma)
    else
        rho, u, p = 0.0, 0.0, 0.0
    end
    
    # reconstruct the fluxesbased on the star state
    rhou = rho * u
    E = p / (eq.gamma - 1.0) + 0.5 * rho * u^2
    
    return [rhou, rhou * u + p, u * (E + p)]
end

function compute_dt(mesh::Mesh1D, eq::Euler1D, values::Matrix{Float64}, CFL::Float64)::Float64
    dx = minimum(mesh.cells_center)
    return CFL * dx / max_wave_speed(eq, values, mesh)
end
