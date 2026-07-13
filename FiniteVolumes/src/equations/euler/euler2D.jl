struct Euler2D{T<:Real, F<:AbstractIdealGasNumericalFlux} <: AbstractEquation2D
    gamma::T
    numerical_flux::F
end

num_vars(::Euler2D) = 4

function flux(eq::Euler2D, U_l::AbstractVector{T}, U_r::AbstractVector{T}, normal::SVector{2, T}) where {T<:Real}

    nx, ny = normal # cosθ, sinθ in Toro p.105

    u_hat_l = nx * U_l[2] + ny * U_l[3]
    v_hat_l =-ny * U_l[2] + nx * U_l[3]

    u_hat_r = nx * U_r[2] + ny * U_r[3]
    v_hat_r =-ny * U_r[2] + nx * U_r[3]

    U_hat_l = SVector(U_l[1], u_hat_l, v_hat_l, U_l[4])
    U_hat_r = SVector(U_r[1], u_hat_r, v_hat_r, U_r[4])

    F = get_flux_2D(eq.numerical_flux, U_hat_l, U_hat_r, eq.gamma)

    return SVector(
                   F[1],
                   nx * F[2] - ny * F[3],
                   ny * F[2] + nx * F[3],
                   F[4]
                  )
end

function compute_dt(mesh::Mesh2D, eq::Euler2D, U::AbstractVector{SVector{4, T}}, CFL::T)::T where {T<:Real}
    _compute_dt_local(u, area, perim) = begin
        rho, u, v, p = _cons_to_prim_euler_ideal_gas(u, eq.gamma)
        (rho<=0 || p<=0) && return T(NaN) # check positivity of pressure and density
        a = sqrt(eq.gamma * p / rho)
        s = sqrt(u*u + v*v) + a
        return 2 * area / (perim * s)
    end
    
    dt_local = _compute_dt_local.(U,mesh.cell_measure, mesh.cell_perimeters)
    return CFL * minimum(dt_local)
end

function output_fields(eq::Euler2D)
    γ = eq.gamma
    [
        OutputField("Rho",  :scalar, U -> U[1]),
        OutputField("U", :vector, U -> SVector(U[2] / U[1], U[3] / U[1])),
        OutputField("p", :scalar, U -> (γ - 1) * (U[4] - (U[2]^2 + U[3]^2) / (2 * U[1]))),
        OutputField("E",   :scalar, U -> U[4]),
        OutputField("e",   :scalar, U -> (U[4] - (U[2]^2 + U[3]^2) / (2 * U[1])) / U[1]),
        OutputField("Momentum", :vector, U -> SVector(U[2], U[3])),
    ]
end

