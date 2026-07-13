"""
helper function to find pressure in the exact riemann solver
"""
function _newton_raphson(f, df, p0; max_it_p::Int = 10, p_tol::T = T(1e-6)) where {T<:Real}
    p = p0
    for _ in 1:max_it_p
        p -= f(p)/df(p) 
        if 2*abs(p-p0)/(p+p0) < p_tol
            return p
        end
        p0 = p
    end
    # error("Pressure didn't converge.")
    # should I do somethgin like if device==:cpu ?
end

"""
Return Star Values for a given Riemann problem
"""
function get_star_values(WL::AbstractVector{T}, WR::AbstractVector{T}, gamma::T; max_it_p::Int=10, p_tol::T=1e-6) where {T<:Real}
    rhoL, uL, pL = WL
    rhoR, uR, pR = WR

    aL = sqrt(gamma*pL/rhoL)
    aR = sqrt(gamma*pR/rhoR)

    # compute p* and u* using Newton
    # build f

    AL = 2 / ((gamma+1) * rhoL)
    BL = ((gamma-1) / (gamma+1)) * pL
    AR = 2 / ((gamma+1) * rhoR)
    BR = ((gamma-1) / (gamma+1)) * pR

    fL(p) = p>pL ? (p-pL)*(sqrt(AL/(p+BL))) : ((2*aL)/(gamma-1)) * ((p/pL)^((gamma-1)/(2*gamma))-1)
    fR(p) = p>pR ? (p-pR)*(sqrt(AR/(p+BR))) : ((2*aR)/(gamma-1)) * ((p/pR)^((gamma-1)/(2*gamma))-1)
    f(p) = fL(p) + fR(p) + (uR - uL)

    # build f'
    dfL(p) = p>pL ? sqrt(AL/(BL+p)) * (1 - (p-pL)/(2*(BL+p))) : 1/(rhoL*aL) * (p/pL) ^ (-(gamma + 1)/(2*gamma))
    dfR(p) = p>pR ? sqrt(AR/(BR+p)) * (1 - (p-pR)/(2*(BR+p))) : 1/(rhoR*aR) * (p/pR) ^ (-(gamma + 1)/(2*gamma))

    df(p) = dfL(p) + dfR(p)

    # choose p0, the starting guess
    # from Toro, we choose the two shock approximation as it seemed to give the best results overall
    pPV = (pL + pR)/2 - (uR - uL) * (rhoL + rhoR) * (aL + aR) / 8
    phat = max(p_tol, pPV)

    gL = sqrt(AL / (phat + BL))
    gR = sqrt(AR / (phat + BR))
    pTS = (gL*pL + gR*pR - (uR-uL)) / (gL + gR)

    p0 = max(p_tol, pTS)
    

    # Newton's method to find p*
    pstar = _newton_raphson(f, df, p0, max_it_p=max_it_p, p_tol=p_tol)

    # compute u*
    ustar = (uL+uR+fR(pstar)-fL(pstar)) / 2

    if pstar > pL
        rhostarL = rhoL * (pstar/pL + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)* pstar/pL + 1)
    else
        rhostarL = rhoL*(pstar/pL)^(1/gamma) 
    end

    if pstar>pR
        rhostarR = rhoR * (pstar/pR + (gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)* pstar/pR + 1)
    else
        rhostarR = rhoR*(pstar/pR)^(1/gamma) 
    end

    return SVector(pstar, ustar, rhostarL, rhostarR)
end

"""
Exact Riemann Solver for the 1D Euler system
"""
function solve_riemann_exact(Xi::T, WL::AbstractVector{T}, WR::AbstractVector{T}, gamma::T; max_it_p::Int=10, p_tol::T=1e-6) where {T<:Real}
    rhoL, uL, pL = WL
    rhoR, uR, pR = WR

    aL = sqrt(gamma*pL/rhoL)
    aR = sqrt(gamma*pR/rhoR)

    pstar, ustar, rhostarL, rhostarR = get_star_values(WL, WR, gamma, max_it_p=max_it_p, p_tol=p_tol)

    # find right case : 
    if Xi < ustar
        # we are on the left of the contact wave
        if pstar>pL     # left shock wave
            SL = uL - aL* sqrt(((gamma+1)/(2*gamma) * pstar/pL + (gamma-1)/(2*gamma)))
            if Xi < SL
                return WL
            else
                return SVector(rhostarL, ustar, pstar)
            end

        else            # left fan
            SHL = uL - aL
            astarL = aL*(pstar/pL)^((gamma-1)/(2*gamma))
            STL = ustar - astarL
            if Xi<SHL
                return WL
            elseif Xi>STL
                return SVector(rhostarL, ustar, pstar)
            else
                rhoLfan = rhoL * ((2/(gamma+1) + (gamma-1)/((gamma+1)*aL) * (uL - Xi))^(2/(gamma-1)) )
                uLfan = 2/(gamma+1) * (aL + (gamma-1)*uL/2 + Xi)
                pLfan = pL * ((2/(gamma+1) + (gamma-1)/((gamma+1)*aL) * (uL - Xi))^(2*gamma/(gamma-1)))
                return SVector(rhoLfan, uLfan, pLfan)
            end
        end
    else
        # we are on the right of the contact wave
        if pstar>pR     # right shock wave
            SR = uR + aR* sqrt(((gamma+1)/(2*gamma) * pstar/pR + (gamma-1)/(2*gamma)))
            if Xi>SR
                return WR
            else
                return SVector(rhostarR, ustar, pstar)
            end
        else            # left fan
            SHR = uR + aR
            astarR = aR*(pstar/pR)^((gamma-1)/(2*gamma))
            STR = ustar + astarR
            if Xi>SHR
                return WR
            elseif Xi<STR
                return SVector(rhostarR, ustar, pstar)
            else
                rhoRfan = rhoR * ((2/(gamma+1) - (gamma-1)/((gamma+1)*aR) * (uR - Xi))^(2/(gamma-1)) )
                uRfan = 2/(gamma+1) * (-aR + (gamma-1)*uR/2 + Xi)
                pRfan = pR * ((2/(gamma+1) - (gamma-1)/((gamma+1)*aR) * (uR - Xi))^(2*gamma/(gamma-1)))
                return SVector(rhoRfan, uRfan, pRfan)
            end
        end
    end
end

