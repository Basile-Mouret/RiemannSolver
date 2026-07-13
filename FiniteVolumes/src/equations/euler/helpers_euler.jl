function _cons_to_prim_euler_ideal_gas(U::SVector{3, T}, gamma::T) where {T<:Real}
    rho = U[1]
    rhou = U[2]
    E = U[3]
    u = rhou / rho
    p = (gamma - 1) * (E - rho * u * u / 2)
    return SVector(rho, u, p)
end

function _cons_to_prim_euler_ideal_gas(U::SVector{4, T}, gamma::T) where {T<:Real}
    rho = U[1]
    rhou = U[2]
    rhov = U[3]
    E = U[4]
    u = rhou / rho
    v = rhov / rho
    p = (gamma - 1) * (E - rho * (u * u + v * v) / 2)
    return SVector(rho, u, v, p)
end

function _prim_to_cons_euler_ideal_gas(W::SVector{3, T}, gamma::T) where {T<:Real}
    rho = W[1]
    u = W[2]
    p = W[3]
    rhou = rho * u
    E = p / (gamma - 1) + rho * u * u / 2
    return SVector(rho, rhou, E)
end

function _prim_to_cons_euler_ideal_gas(W::SVector{4, T}, gamma::T) where {T<:Real}
    rho = W[1]
    u = W[2]
    v = W[3]
    p = W[4]
    rhou = rho * u
    rhov = rho * v
    E = p / (gamma - 1) + rho * (u * u + v * v) / 2
    return SVector(rho, rhou, rhov, E)
end

