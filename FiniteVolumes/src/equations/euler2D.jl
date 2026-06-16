struct Euler2D <: AbstractEquation2D
    gamma::Float64
    solver::Symbol
end

function _cons_to_prim_euler_2D_ideal_gas(U::AbstractVector{T}, gamma::T) where {T<:Real}
    rho = U[1]
    rhou = U[2]
    rhov = U[3]
    E = U[4]
    u = rhou / rho
    v = rhov / rho
    p = (gamma - 1.0) * (E - 0.5 * (rhou^2+rhov^2) / rho)
    return SVector(rho, u, v, p)
end

function _prim_to_cons_euler_2D_ideal_gas(W::AbstractVector{T}, gamma::T) where {T<:Real}
    rho = W[1]
    u = W[2]
    v = W[3]
    p = W[4]
    rhou = rho * u
    rhov = rho * v
    E = p / (gamma - 1.0) + 0.5 * rho * (u^2+v^2)
    return SVector(rho, rhou, rhov, E)
end

num_vars(::Euler2D) = 4

function max_wave_speed(::Mesh2D, eq::Euler2D, U::Matrix{Float64})
    s_max = 0.0
    nvars = num_vars(eq)
    for i in axes(U, 1)
        U_cell = SVector{nvars}(@view U[i, :])

        rho, u, v, p = _cons_to_prim_euler_2D_ideal_gas(U_cell, eq.gamma)

        a = sqrt(eq.gamma * p / rho)
        s_max = max(s_max, sqrt(u*u + v*v) + a)
    end

    return s_max
end


function flux(eq::Euler2D, UL::AbstractVector{Float64}, UR::AbstractVector{Float64}, normal::NTuple{2, Float64})

    nx, ny = normal # cosθ, sinθ in Toro p.105
    # tangantial component

    WL = _cons_to_prim_euler_2D_ideal_gas(UL, eq.gamma)
    WR = _cons_to_prim_euler_2D_ideal_gas(UR, eq.gamma)


    # solve the Riemann problem at the interface (x/t = 0.0)
    if eq.solver == :exact
        rho, uh, p = solve_riemann_exact(0.0,
                                         SVector(WL[1], WL[2] * nx + WL[3] * ny, WL[4]), #projecting WL and WR on the normal
                                         SVector(WR[1], WR[2] * nx + WR[3] * ny, WR[4]),
                                         eq.gamma)
    else
        error("This Rieman solver is not implemented")
    end

    #tangential component is chosen based on the upwind scheme
    vh = uh>=0.0 ? WL[2] * (-ny) + WL[3] * nx : WR[2] * (-ny) + WR[3] * nx

    # reconstruct the fluxes based on the star state
    E = p / (eq.gamma - 1.0) + 0.5 * rho * (uh^2 + vh^2)

    return SVector(rho*uh,
            nx*(rho*uh^2 + p) - ny*(rho*uh*vh),
            ny*(rho*uh^2 + p) + nx*(rho*uh*vh),
            uh*(E+p))
end

function compute_dt(mesh::Mesh2D, eq::Euler2D, values::Matrix{Float64}, CFL::Float64)::Float64
    return CFL * 2.0 * minimum(mesh.cell_measure[i] / mesh.cell_perimeters[i] for i in eachindex(mesh.cell_measure)) / max_wave_speed(mesh, eq, values)
end

