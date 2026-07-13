struct IdealGasHLL <: AbstractIdealGasNumericalFlux
end

"""
Compute the HLL flux in 3D
"""
function get_flux_3D(numerical_flux::IdealGasHLL, U_l::AbstractVector{T}, U_r::AbstractVector{T}, gamma::T) where {T<:Real}
    # compute the wave speeds S_l and S_r
    # using Einfeldt estimates based on Roe eigenvalues
    rho_l, rhou_l, rhov_l, rhow_l, E_l = U_l
    u_l = rhou_l / rho_l
    v_l = rhov_l / rho_l
    w_l = rhow_l / rho_l
    p_l = (gamma-1) * (E_l - 0.5 * rho_l * (u_l*u_l + v_l*v_l + w_l*w_l))
    a_l = sqrt(gamma * p_l / rho_l)
    sqrt_rho_l = sqrt(rho_l)

    rho_r, rhou_r, rhov_r, rhow_r, E_r = U_r
    u_r = rhou_r / rho_r
    v_r = rhov_r / rho_r
    w_r = rhow_r / rho_r
    p_r = (gamma-1) * (E_r - 0.5 * rho_r * (u_r*u_r + v_r*v_r + w_r*w_r))
    a_r = sqrt(gamma * p_r / rho_r)
    sqrt_rho_r = sqrt(rho_r)

    d_bar = sqrt((sqrt_rho_l*a_l*a_l + sqrt_rho_r*a_r*a_r)/(sqrt_rho_l + sqrt_rho_r) + (0.5  * (sqrt_rho_l * sqrt_rho_r)/((sqrt_rho_l + sqrt_rho_r)^2)) * (u_r -u_l)^2)

    u_bar =  (sqrt_rho_l * u_l + sqrt_rho_r * u_r) / (sqrt_rho_l + sqrt_rho_r)

    S_l = min(u_l - a_l, u_bar - d_bar)
    S_r = max(u_r + a_r, u_bar + d_bar)

    # compute the HLL flux
    
    F_l = SVector(rhou_l, rhou_l*u_l + p_l, rhou_l*v_l, rhou_l*w_l, u_l*(E_l + p_l))
    F_r = SVector(rhou_r, rhou_r*u_r + p_r, rhou_r*v_r, rhou_r*w_r, u_r*(E_r + p_r))

    if S_l >= 0 
        return F_l
    elseif S_r <= 0
        return F_r
    else
        return (S_r*F_l - S_l*F_r + S_l*S_r*(U_r - U_l))/(S_r - S_l)
    end
end

"""
HLL flux for Godunov's scheme in 1D
"""
function get_flux_1D(numerical_flux::IdealGasHLL, U_l::AbstractVector{T}, U_r::AbstractVector{T}, gamma::T) where {T<:Real}
    U_pad_l = SVector(U_l[1], U_l[2], 0, 0, U_l[3])
    U_pad_r = SVector(U_r[1], U_r[2], 0, 0, U_r[3])
    F = get_flux_3D(numerical_flux, U_pad_l, U_pad_r, gamma)
    return SVector(F[1], F[2], F[5])
end

"""
HLL flux for Godunov's scheme in 2D
"""
function get_flux_2D(numerical_flux::IdealGasHLL, U_l::AbstractVector{T}, U_r::AbstractVector{T}, gamma::T) where {T<:Real}
    U_pad_l = SVector(U_l[1], U_l[2], U_l[3], 0, U_l[4])
    U_pad_r = SVector(U_r[1], U_r[2], U_r[3], 0, U_r[4])
    F = get_flux_3D(numerical_flux, U_pad_l, U_pad_r, gamma)
    return SVector(F[1], F[2], F[3], F[5])
end

