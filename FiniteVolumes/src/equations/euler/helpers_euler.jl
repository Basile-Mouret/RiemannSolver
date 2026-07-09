
function _cons_to_prim_euler_ideal_gas(U::AbstractVector{T}, gamma::T, dim::Symbol) where {T<:Real}
    if dim == :one
        rho = U[1]
        rhou = U[2]
        E = U[3]
        u = rhou / rho
        p = (gamma - 1.0) * (E - 0.5 * rho * u * u)
        return SVector(rho, u, p)
    elseif dim == :two
        rho = U[1]
        rhou = U[2]
        rhov = U[3]
        E = U[4]
        u = rhou / rho
        v = rhov / rho
        p = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v))
        return SVector(rho, u, v, p)
    else
        error("Not implemented for this many dimensions")
    end
end

function _prim_to_cons_euler_ideal_gas(W::AbstractVector{T}, gamma::T, dim::Symbol) where {T<:Real}
    if dim == :one
        rho = W[1]
        u = W[2]
        p = W[3]
        rhou = rho * u
        E = p / (gamma - 1.0) + 0.5 * rho * u * u
        return SVector(rho, rhou, E)
    elseif dim == :two
        rho = W[1]
        u = W[2]
        v = W[3]
        p = W[4]
        rhou = rho * u
        rhov = rho * v
        E = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v)
        return SVector(rho, rhou, rhov, E)
    else
        error("Not implemented for this many dimensions")
    end
end

