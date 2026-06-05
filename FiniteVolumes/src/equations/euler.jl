
struct Euler1D <: AbstractEquation1D
    # we could define a solver structure such that there is an exact/approximate solver
end

function _cons_to_prim_euler_1D_ideal_gas(U::Vector{T}, gamma::T) where {T<:Real}
    return [U[1], U[2]/U[1], (gamma-1)*U[1]*(U[3]-0.5*U[2]*U[2]/U[1])]
end

function _prim_to_cons_euler_1D_ideal_gas(W::Vector{T}, gamma::T) where {T<:Real}
    return [W[1], W[2]*W[1], W[1]*(0.5*W[2]*W[2] + 1/((gamma-1)*W[1]) * W[3])]
end

num_vars(::Euler1D) = 4

function max_wave_speed(eq::Euler1D, values::Matrix{Float64}, Mesh::Mesh1D )
    # TODO
end

function flux(eq::Euler1D, UL::Vector{Float64}, UR::Vector{Float64})
    gamma = Nothing
    WL = _cons_to_prim_euler_1D_ideal_gas(UL, gamma)
    WR = _cons_to_prim_euler_1D_ideal_gas(UR, gamma)
    rho, u, p = solve_riemann(0.0, WL, WR, gamma, aL, aR)
    E = rho*(0.5*u*u + p/((gamma-1)*rho))
    return [rho*u, rho*u*u + p, u*(E+p)]
end

function apply_ghost(bc::Reflecting, u_interior::Vector{Float64}, ::Float64)
    # TODO
end

