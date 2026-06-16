struct Advection1D <: AbstractEquation1D
    c::Float64
    flux_type::Symbol
end

Advection1D(; c::Float64, flux_type::Symbol = :upwind) = Advection1D(c, flux_type)

num_vars(::Advection1D) = 1

function flux(eq::Advection1D, uL::AbstractVector{Float64}, uR::AbstractVector{Float64})
    if eq.flux_type == :upwind
        F = eq.c > 0 ? eq.c * uL[1] : eq.c * uR[1]
    elseif eq.flux_type == :centered
        F = 0.5 * eq.c * (uL[1] + uR[1])
    elseif eq.flux_type == :left
        F = eq.c * uL[1]
    elseif eq.flux_type == :right
        F = eq.c * uR[1]
    else
        error("Unknown flux_type: $(eq.flux_type)")
    end
    return SVector(F)
end

function exact_solution!(utrue::Matrix{Float64}, eq::Advection1D, xmid::Vector{Float64},
                         ic::Function, bcs::Dict, x0::Float64, x1::Float64, t::Float64)
    c = eq.c
    if bcs[:left] isa Outflow && bcs[:right] isa Outflow
        for i in eachindex(xmid)
            utrue[i, 1] = ic(mod(xmid[i] - c * t, x1 - x0))[1]
        end
    else
        if c > 0
            for i in eachindex(xmid)
                x = xmid[i]
                if x <= c * t
                    arrival_time = t - x / c
                    bcs[:left] isa Dirichlet || error("inflow BC must be Dirichlet")
                    utrue[i, 1] = bcs[:left].value(arrival_time)
                else
                    utrue[i, 1] = ic(x - c * t)[1]
                end
            end
        else
            for i in eachindex(xmid)
                x = xmid[i]
                if x >= x1 + c * t
                    arrival_time = t + (x1 - x) / c
                    bcs[:right] isa Dirichlet || error("inflow BC must be Dirichlet")
                    utrue[i, 1] = bcs[:right].value(arrival_time)
                else
                    utrue[i, 1] = ic(x - c * t)[1]
                end
            end
        end
    end
end

function compute_dt(mesh::Mesh1D, eq::Advection1D, values::Matrix{Float64}, CFL::Float64)::Float64
    dx = minimum(mesh.cell_measure)
    return CFL * dx / abs(eq.c)
end

