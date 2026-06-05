struct Euler1D <: AbstractEquation1D
    max_wave_speed::Float64
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
    gamma = 5.0 / 3.0
    s_max = 0.0
    
    # Iterate over all cells to find the fastest global wave speed
    # Assuming U is a matrix where each column U[:, i] is the state of cell i
    for i in axes(U, 2)
        U_cell = U[:, i]
        
        # Convert to primitive to easily get u and p
        rho, u, p = _cons_to_prim_euler_1D_ideal_gas(U_cell, gamma)
        
        # Calculate local speed of sound
        a = sqrt(gamma * p / rho)
        
        # Update the global maximum wave speed
        s_max = max(s_max, abs(u) + a)
    end
    
    return s_max
end

function flux(::Euler1D, UL::Vector{Float64}, UR::Vector{Float64})
    gamma = 5.0 / 3.0         # gaz monoatomique
    
    WL = _cons_to_prim_euler_1D_ideal_gas(UL, gamma)
    WR = _cons_to_prim_euler_1D_ideal_gas(UR, gamma)
    
    # solve the Riemann problem at the interface (x/t = 0.0)
    rho, u, p = solve_riemann(0.0, WL, WR, gamma)
    
    # reconstruct the fluxesbased on the star state
    rhou = rho * u
    E = p / (gamma - 1.0) + 0.5 * rho * u^2
    
    return [rhou, rhou * u + p, u * (E + p)]
end

