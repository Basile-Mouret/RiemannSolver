"""
helper function to find pressure in the exact riemann solver
"""
function _newton_raphson(f, df, p0; max_it::Int = 10, TOL::Float64 = 10^(-6))
    p = p0
    for _ in 1:max_it
        p -= f(p)/df(p) 
        if 2.0*abs(p-p0)/(p+p0) < TOL
            return p
        end
        p0 = p
    end
    return p
end

"""
Exact Riemann Solver for the 1D Euler system
"""
function solve_riemann_exact(S::T, WL::AbstractVector{T}, WR::AbstractVector{T}, gamma::T, ; TOL::T=1e-6) where {T<:Real}
    rhoL, uL, pL = WL
    rhoR, uR, pR = WR

    aL = sqrt(gamma*pL/rhoL)
    aR = sqrt(gamma*pR/rhoR)

    # compute p* and u* using Newton
    # build f

    AL = 2.0 / ((gamma+1.0) * rhoL)
    BL = ((gamma-1.0) / (gamma+1.0)) * pL
    AR = 2.0 / ((gamma+1.0) * rhoR)
    BR = ((gamma-1.0) / (gamma+1.0)) * pR

    fL(p) = p>pL ? (p-pL)*(sqrt(AL/(p+BL))) : ((2.0*aL)/(gamma-1.0)) * ((p/pL)^((gamma-1.0)/(2.0*gamma))-1.0)
    fR(p) = p>pR ? (p-pR)*(sqrt(AR/(p+BR))) : ((2.0*aR)/(gamma-1.0)) * ((p/pR)^((gamma-1.0)/(2.0*gamma))-1.0)
    f(p) = fL(p) + fR(p) + (uR - uL)

    # build f'
    dfL(p) = p>pL ? sqrt(AL/(BL+p)) * (1.0 - (p-pL)/(2.0*(BL+p))) : 1.0/(rhoL*aL) * (p/pL) ^ (-(gamma + 1.0)/(2.0*gamma))
    dfR(p) = p>pR ? sqrt(AR/(BR+p)) * (1.0 - (p-pR)/(2.0*(BR+p))) : 1.0/(rhoR*aR) * (p/pR) ^ (-(gamma + 1.0)/(2.0*gamma))

    df(p) = dfL(p) + dfR(p)

    # choose p0, the starting guess
    # from Toro, we choose the two shock approximation as it seemed to give the best results overall
    pPV = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (rhoL + rhoR) * (aL + aR)
    phat = max(TOL, pPV)

    gL = sqrt(AL / (phat + BL))
    gR = sqrt(AR / (phat + BR))
    pTS = (gL*pL + gR*pR - (uR-uL)) / (gL + gR)

    p0 = max(TOL, pTS)
    

    # Newton's method to find p*
    pstar = _newton_raphson(f, df, p0)

    # compute u*
    ustar = 0.5*(uL+uR+fR(pstar)-fL(pstar)) 
    


    # find right case : 
    if S < ustar
        # we are on the left of the contact wave
        if pstar>pL     # left shock wave
            SL = uL - aL* sqrt(((gamma+1.0)/(2*gamma) * pstar/pL + (gamma-1.0)/(2.0*gamma)))
            if S < SL
                return WL
            else
                rhostarL = rhoL * (pstar/pL + (gamma-1.0)/(gamma+1.0))/((gamma-1.0)/(gamma+1.0)* pstar/pL + 1.0)
                return SVector(rhostarL, ustar, pstar)
            end

        else            # left fan
            SHL = uL - aL
            astarL = aL*(pstar/pL)^((gamma-1.0)/(2.0*gamma))
            STL = ustar - astarL
            if S<SHL
                return WL
            elseif S>STL
                rhostarL = rhoL*(pstar/pL)^(1.0/gamma) 
                return SVector(rhostarL, ustar, pstar)
            else
                rhoLfan = rhoL * ((2.0/(gamma+1.0) + (gamma-1.0)/((gamma+1.0)*aL) * (uL - S))^(2.0/(gamma-1.0)) )
                uLfan = 2.0/(gamma+1.0) * (aL + 0.5*(gamma-1.0)*uL + S)
                pLfan = pL * ((2.0/(gamma+1.0) + (gamma-1.0)/((gamma+1.0)*aL) * (uL - S))^(2.0*gamma/(gamma-1.0)))
                return SVector(rhoLfan, uLfan, pLfan)
            end
        end
    else
        # we are on the right of the contact wave
        if pstar>pR     # right shock wave
            SR = uR + aR* sqrt(((gamma+1.0)/(2*gamma) * pstar/pR + (gamma-1.0)/(2.0*gamma)))
            if S>SR
                return WR
            else
                rhostarR = rhoR * (pstar/pR + (gamma-1.0)/(gamma+1.0))/((gamma-1.0)/(gamma+1.0)* pstar/pR + 1.0)
                return SVector(rhostarR, ustar, pstar)
            end
        else            # left fan
            SHR = uR + aR
            astarR = aR*(pstar/pR)^((gamma-1.0)/(2.0*gamma))
            STR = ustar + astarR
            if S>SHR
                return WR
            elseif S<STR
                rhostarR = rhoR*(pstar/pR)^(1.0/gamma) 
                return SVector(rhostarR, ustar, pstar)
            else
                rhoRfan = rhoR * ((2.0/(gamma+1.0) - (gamma-1.0)/((gamma+1.0)*aR) * (uR - S))^(2.0/(gamma-1.0)) )
                uRfan = 2.0/(gamma+1.0) * (-aR + 0.5*(gamma-1.0)*uR + S)
                pRfan = pR * ((2.0/(gamma+1.0) - (gamma-1.0)/((gamma+1.0)*aR) * (uR - S))^(2.0*gamma/(gamma-1.0)))
                return SVector(rhoRfan, uRfan, pRfan)
            end
        end
    end
end

"""
Roe
"""
function solve_riemann_roe(S::T, WL::Vector{T}, WR::Vector{T}, gamma::T, ; TOL::T=1e-6) where {T<:Real}
end

"""
HLL
"""
function solve_riemann_hll(S::T, WL::Vector{T}, WR::Vector{T}, gamma::T, ; TOL::T=1e-6) where {T<:Real}
    
end

"""
HLLC
"""
function solve_riemann_hllc(S::T, WL::Vector{T}, WR::Vector{T}, gamma::T, ; TOL::T=1e-6) where {T<:Real}
end
