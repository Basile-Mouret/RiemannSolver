struct Euler1D <: AbstractEquation1D
    gamma::Float64
    numerical_flux::Symbol
end

num_vars(::Euler1D) = 3

function flux(eq::Euler1D, UL::AbstractVector{Float64}, UR::AbstractVector{Float64})
    if eq.numerical_flux == :Godunov
        return get_godunov_flux_1D(UL, UR, eq.gamma)
    elseif eq.numerical_flux == :Roe
        return get_roe_flux_1D(UL, UR, eq.gamma)
    else
        error("flux not implemented")
    end
end

function compute_dt(mesh::Mesh1D, eq::Euler1D, U::Matrix{Float64}, CFL::Float64)::Float64
    dt_min = Inf
    nvars = num_vars(eq)

    for i in axes(U, 1)
        U_cell = SVector{nvars}(@view U[i, :])
        rho, u, p = _cons_to_prim_euler_ideal_gas(U_cell, eq.gamma, :one)
        a = sqrt(eq.gamma * p / rho)
        dx = mesh.cell_measure[i]
        dt_min = min(dt_min, dx/(abs(u) + a))
    end

    return CFL * dt_min
end

function output_fields(eq::Euler1D)
    γ = eq.gamma
    [
        OutputField("Rho", :scalar, U -> U[1]),
        OutputField("U", :scalar, U -> U[2] / U[1]),
        OutputField("p", :scalar, U -> (γ - 1.0) * (U[3] - 0.5 * U[2]^2 / U[1])),
        OutputField("E", :scalar, U -> U[3]),
        OutputField("e", :scalar, U -> (U[3] - 0.5 * U[2]^2 / U[1])/ U[1]),
    ]
end

