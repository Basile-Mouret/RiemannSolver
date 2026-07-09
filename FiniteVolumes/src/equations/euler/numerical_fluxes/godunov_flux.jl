struct IdealGasGodunov <: AbstractIdealGasNumericalFlux
    # initial_p_guess::Symbol
    # max_it_newton::Int
    
end

function get_flux_1D(numerical_flux::IdealGasGodunov, U_l::AbstractVector{T}, U_r::AbstractVector{T}, gamma::T) where {T<:Real}
        # get the primitive variables
        W_l = _cons_to_prim_euler_ideal_gas(U_l, gamma, :one)
        W_r = _cons_to_prim_euler_ideal_gas(U_r, gamma, :one)

        # use an exact riemann solver
        rho, uh, p = solve_riemann_exact(0.0, W_l, W_r, gamma)

        E = p / (gamma - 1.0) + 0.5 * rho * (uh^2)

        # return numerical flux
        return SVector(rho*uh, rho*uh^2 + p, uh*(E+p))
end

function get_flux_2D(numerical_flux::IdealGasGodunov,U_l::AbstractVector{T}, U_r::AbstractVector{T}, gamma::T) where {T<:Real}
        # get the primitive variables
        W_l = _cons_to_prim_euler_ideal_gas(U_l, gamma, :two)
        W_r = _cons_to_prim_euler_ideal_gas(U_r, gamma, :two)

        # use an exact riemann solver
        rho, u_h, p = solve_riemann_exact(0.0,
                                          SVector(W_l[1], W_l[2], W_l[4]),
                                          SVector(W_r[1], W_r[2], W_r[4]),
                                          gamma)
        # get tangent component from the orientation of the contact wave
        v_h = u_h>0 ? W_l[3] : W_r[3]
        E = p / (gamma - 1.0) + 0.5*rho*(u_h*u_h + v_h*v_h)
        return SVector(rho*u_h, rho*u_h^2 + p, rho*u_h*v_h, u_h*(E+p))
end

