struct Euler2D <: AbstractEquation2D
    gamma::Float64
    numerical_flux::Symbol
end

num_vars(::Euler2D) = 4

function flux(eq::Euler2D, U_l::AbstractVector{Float64}, U_r::AbstractVector{Float64}, normal::SVector{2, Float64})

    nx, ny = normal # cosθ, sinθ in Toro p.105

    u_hat_l = nx * U_l[2] + ny * U_l[3]
    v_hat_l =-ny * U_l[2] + nx * U_l[3]

    u_hat_r = nx * U_r[2] + ny * U_r[3]
    v_hat_r =-ny * U_r[2] + nx * U_r[3]

    U_hat_l = SVector(U_l[1], u_hat_l, v_hat_l, U_l[4])
    U_hat_r = SVector(U_r[1], u_hat_r, v_hat_r, U_r[4])

    if eq.numerical_flux == :Godunov
        F = get_godunov_flux_2D(U_hat_l, U_hat_r, eq.gamma)
    elseif eq.numerical_flux == :HLL
        F = get_hllc_flux_2D(U_hat_l, U_hat_r, eq.gamma)
    elseif eq.numerical_flux == :HLLC
        F = get_hllc_flux_2D(U_hat_l, U_hat_r, eq.gamma)
    elseif eq.numerical_flux == :Roe
        F = get_roe_flux_2D(U_hat_l, U_hat_r, eq.gamma)
    else
        error("This solver is not implemented")
    end

    return SVector(
                   F[1],
                   nx * F[2] - ny * F[3],
                   ny * F[2] + nx * F[3],
                   F[4]
                  )
end

function compute_dt(mesh::Mesh2D, eq::Euler2D, U::Matrix{Float64}, CFL::Float64)::Float64
    # not optimal, I could try to efficiently compute minimum(cell_meas/(s_max_into_cell * cell_perimeter))
    # this way I can have finer mesh at low speed locations

    dt_min = Inf
    nvars = num_vars(eq)

    for i in axes(U, 1)
        U_cell = SVector{nvars}(@view U[i, :])
        rho, u, v, p = _cons_to_prim_euler_ideal_gas(U_cell, eq.gamma, :two)

        a = sqrt(eq.gamma * p / rho) 
        s_local = sqrt(u*u + v*v) + a

        char_length = 2.0 * mesh.cell_measure[i] / mesh.cell_perimeters[i]
        dt_local = char_length / s_local
        dt_min = min(dt_min, dt_local)
    end

    return CFL * dt_min
end

function output_fields(eq::Euler2D)
    γ = eq.gamma
    [
        OutputField("Rho",  :scalar, U -> U[1]),
        OutputField("U", :vector, U -> SVector(U[2] / U[1], U[3] / U[1])),
        OutputField("p", :scalar, U -> (γ - 1.0) * (U[4] - 0.5 * (U[2]^2 + U[3]^2) / U[1])),
        OutputField("E",   :scalar, U -> U[4]),
        OutputField("e",   :scalar, U -> (U[4] - 0.5 * (U[2]^2 + U[3]^2) / U[1]) / U[1]),
        OutputField("Momentum", :vector, U -> SVector(U[2], U[3])),
    ]
end

