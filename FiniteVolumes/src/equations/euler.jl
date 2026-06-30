struct Euler1D <: AbstractEquation1D
    gamma::Float64
    solver::Symbol
end

function _cons_to_prim_euler_1D_ideal_gas(U::AbstractVector{T}, gamma::T) where {T<:Real}
    rho = U[1]
    rhou = U[2]
    E = U[3]
    u = rhou / rho
    p = (gamma - 1.0) * (E - 0.5 * rhou^2 / rho)
    return SVector(rho, u, p)
end

function _prim_to_cons_euler_1D_ideal_gas(W::AbstractVector{T}, gamma::T) where {T<:Real}
    rho = W[1]
    u = W[2]
    p = W[3]
    rhou = rho * u
    E = p / (gamma - 1.0) + 0.5 * rho * u^2
    return SVector(rho, rhou, E)
end

num_vars(::Euler1D) = 3

function max_wave_speed(::Mesh1D, eq::Euler1D, U::Matrix{Float64})
    s_max = 0.0
    nvars = num_vars(eq)

    for i in axes(U, 1)
        U_cell = SVector{nvars}(@view U[i, :])

        rho, u, p = _cons_to_prim_euler_1D_ideal_gas(U_cell, eq.gamma)

        a = sqrt(eq.gamma * p / rho)
        s_max = max(s_max, abs(u) + a)
    end

    return s_max
end

function flux(eq::Euler1D, UL::AbstractVector{Float64}, UR::AbstractVector{Float64})

    WL = _cons_to_prim_euler_1D_ideal_gas(UL, eq.gamma)
    WR = _cons_to_prim_euler_1D_ideal_gas(UR, eq.gamma)

    function F(W)
        rho, u, p = W
        rhou = rho * u
        E = p / (eq.gamma - 1.0) + 0.5 * rho * u^2
        return SVector(rhou, rhou * u + p, u * (E + p))
    end

    # solve the Riemann problem at the interface (x/t = 0.0)
    if eq.solver == :exact
        Ws = solve_riemann_exact(0.0, WL, WR, eq.gamma)
    # elseif eq.solver == :Lax_Friedrich
    #     return 0.5*(F(WL) + F(WR)) + 0.5*(dt/dx)*(UL-UR) I need dt and dx in my flux call
    elseif eq.solver == :hll
        
    else
        rho, u, p = 0.0, 0.0, 0.0
    end

    # reconstruct the fluxesbased on the star state

    return F(Ws) 
end

function compute_dt(mesh::Mesh1D, eq::Euler1D, values::Matrix{Float64}, CFL::Float64)::Float64
    dx = minimum(mesh.cell_measure)
    return CFL * dx / max_wave_speed(mesh, eq, values)
end

function output_fields(eq::Euler1D)
    γ = eq.gamma
    [
        OutputField("density",  :scalar, U -> U[1]),
        OutputField("velocity", :scalar, U -> U[2] / U[1]),
        OutputField("pressure", :scalar, U -> (γ - 1.0) * (U[3] - 0.5 * U[2]^2 / U[1])),
        OutputField("Energy",   :scalar, U -> U[3]),
    ]
end


## boundary conditions


