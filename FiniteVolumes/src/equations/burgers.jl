struct Burgers1D <: AbstractEquation1D end

num_vars(::Burgers1D) = 1

function flux(::Burgers1D, Ul::Vector{Float64}, Ur::Vector{Float64})
    f(x) = 0.5*x*x
    ul, ur = Ul[1], Ur[1]
    if ul > ur
        S = 0.5*(ul+ur)
        if S> 0
            return f(ul)
        else 
            return f(ur)
        end
    else
        if ul >= 0
            return f(ul)
        elseif ur <= 0
            return f(ur)
        else
            return f(0.)
        end
    end
end

function compute_dt(mesh::Mesh1D, eq::Burgers1D, values::Matrix{Float64}, CFL::Float64)::Float64
    return minimum(CFL .* mesh.cells_length ./ abs.(values))
end
