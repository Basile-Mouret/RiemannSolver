"""
Original Roe flux for a 1D Riemann Problem
"""
function get_roe_flux_1D(U_l::AbstractVector{T}, U_r::AbstractVector{T}, gamma::T) where {T<:Real}
    U_pad_l = SVector(U_l[1], U_l[2], 0, 0, U_l[3])
    U_pad_r = SVector(U_r[1], U_r[2], 0, 0, U_r[3])
    F = get_roe_flux_3D(U_pad_l, U_pad_r, gamma)
    return SVector(F[1], F[2], F[5])
end

"""
Original Roe flux for a 2D Riemann Problem
"""
function get_roe_flux_2D(U_l::AbstractVector{T}, U_r::AbstractVector{T}, gamma::T) where {T<:Real}
    U_pad_l = SVector(U_l[1], U_l[2], U_l[3], 0, U_l[4])
    U_pad_r = SVector(U_r[1], U_r[2], U_r[3], 0, U_r[4])
    F = get_roe_flux_3D(U_pad_l, U_pad_r, gamma)
    return SVector(F[1], F[2], F[3], F[5])
end

"""
Original Roe flux for a 3D Riemann Problem
"""
function get_roe_flux_3D(U_l::AbstractVector{T}, U_r::AbstractVector{T}, gamma::T) where {T<:Real}
    rho_l, rhou_l, rhov_l, rhow_l, E_l = U_l
    u_l = rhou_l / rho_l
    v_l = rhov_l / rho_l
    w_l = rhow_l / rho_l
    p_l = (gamma-1) * (E_l - 0.5 * rho_l * (u_l*u_l + v_l*v_l + w_l*w_l))
    H_l = (E_l + p_l) / rho_l

    rho_r, rhou_r, rhov_r, rhow_r, E_r = U_r
    u_r = rhou_r / rho_r
    v_r = rhov_r / rho_r
    w_r = rhow_r / rho_r
    p_r = (gamma-1) * (E_r - 0.5 * rho_r * (u_r*u_r + v_r*v_r + w_r*w_r))
    H_r = (E_r + p_r) / rho_r

    sqrt_rho_l = sqrt(rho_l)
    sqrt_rho_r = sqrt(rho_r)

    # compute the roe average values
    u_tilde =  (sqrt_rho_l * u_l + sqrt_rho_r * u_r) / (sqrt_rho_l + sqrt_rho_r)
    v_tilde =  (sqrt_rho_l * v_l + sqrt_rho_r * v_r) / (sqrt_rho_l + sqrt_rho_r)
    w_tilde =  (sqrt_rho_l * w_l + sqrt_rho_r * w_r) / (sqrt_rho_l + sqrt_rho_r)
    H_tilde =  (sqrt_rho_l * H_l + sqrt_rho_r * H_r) / (sqrt_rho_l + sqrt_rho_r)
    a_tilde = sqrt((gamma-1.0) * (H_tilde - 0.5 * (u_tilde*u_tilde + v_tilde*v_tilde + w_tilde*w_tilde)))

    # compute the averaged eigenvalues
    lam_1 = u_tilde - a_tilde
    lam_2 = u_tilde
    lam_3 = u_tilde
    lam_4 = u_tilde
    lam_5 = u_tilde + a_tilde

    # compute the averaged right eigenvectors
    K_tilde_1 = SVector(1.0, u_tilde - a_tilde, v_tilde, w_tilde, H_tilde - u_tilde * a_tilde)
    K_tilde_2 = SVector(1, u_tilde, v_tilde, w_tilde, 0.5 * (u_tilde*u_tilde + v_tilde*v_tilde + w_tilde*w_tilde))
    K_tilde_3 = SVector(0.0, 0.0, 1.0, 0.0, v_tilde)
    K_tilde_4 = SVector(0.0, 0.0, 0.0, 1.0, w_tilde)
    K_tilde_5 = SVector(1.0, u_tilde + a_tilde, v_tilde, w_tilde, H_tilde + u_tilde * a_tilde)

    # compute the wave strengths
    ΔU = U_r - U_l

    alpha_tilde_3 = (ΔU[3]) - v_tilde*(ΔU[1])
    alpha_tilde_4 = (ΔU[4]) - w_tilde*(ΔU[1])
    alpha_tilde_2 = ((gamma-1) / (a_tilde*a_tilde)) * (
                        (ΔU[1]) * (H_tilde - u_tilde*u_tilde)
                        + u_tilde*(ΔU[2])
                        - (ΔU[5] - alpha_tilde_3*v_tilde - alpha_tilde_4*w_tilde)
                       )
    alpha_tilde_1 = (0.5/a_tilde) * (ΔU[1] * (u_tilde + a_tilde) - ΔU[2] - a_tilde * alpha_tilde_2)
    alpha_tilde_5 = ΔU[1] - (alpha_tilde_1 + alpha_tilde_2)

    # compute the intercell flux
    
    F_l = SVector(rhou_l, rhou_l*u_l + p_l, rhou_l*v_l, rhou_l*w_l, u_l*(E_l + p_l))
    F_r = SVector(rhou_r, rhou_r*u_r + p_r, rhou_r*v_r, rhou_r*w_r, u_r*(E_r + p_r))

    return 0.5*(F_l + F_r) - 0.5*(alpha_tilde_1*abs(lam_1)*K_tilde_1
                                        + alpha_tilde_2*abs(lam_2)*K_tilde_2
                                        + alpha_tilde_3*abs(lam_3)*K_tilde_3
                                        + alpha_tilde_4*abs(lam_4)*K_tilde_4
                                        + alpha_tilde_5*abs(lam_5)*K_tilde_5
                                       )
end
