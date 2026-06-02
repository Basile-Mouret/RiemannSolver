struct Wave1D <: AbstractEquation1D
    kappa::Float64
    rho::Float64
end

num_vars(::Wave1D) = 2

max_wave_speed(eq::Wave1D) = sqrt(eq.kappa / eq.rho)

function flux(eq::Wave1D, uL::Vector{Float64}, uR::Vector{Float64})
    pl, ul = uL[1], uL[2]
    pr, ur = uR[1], uR[2]
    c = sqrt(eq.kappa / eq.rho)
    Fp = (0.5 / eq.rho) * (ur + ul) - 0.5 * c * (pr - pl)
    Fu = 0.5 * eq.kappa * (pr + pl) - 0.5 * c * (ur - ul)
    return [Fp, Fu]
end

function apply_ghost(bc::Reflecting, u_interior::Vector{Float64}, ::Float64)
    return [u_interior[1], -u_interior[2]]
end

function exact_solution!(utrue::Matrix{Float64}, eq::Wave1D, xmid::Vector{Float64},
                         ic::Function, bcs::Dict, x0::Float64, x1::Float64, t::Float64)
    N = length(xmid)
    c = sqrt(eq.kappa / eq.rho)
    rho_c = eq.rho * c
    L = x1 - x0
    if bcs[:left] isa Outflow && bcs[:right] isa Outflow
        for i in 1:N
            pv, uv = ic(mod(xmid[i] + c * t, L))
            v = 0.5 * pv - 0.5 / rho_c * uv
            pw, uw = ic(mod(xmid[i] - c * t, L))
            w = 0.5 * pw + 0.5 / rho_c * uw
            utrue[i, 1] = v + w
            utrue[i, 2] = rho_c * (w - v)
        end
    else
        utrue .= 0.0
    end
end

function entropy(eq::Wave1D, cell_values::Matrix{Float64}, dx::Float64)
    p_sq = sum(cell_values[:, 1] .^ 2)
    u_sq = sum(cell_values[:, 2] .^ 2)
    return 0.5 * dx * (p_sq + u_sq / (eq.rho * eq.kappa))
end
