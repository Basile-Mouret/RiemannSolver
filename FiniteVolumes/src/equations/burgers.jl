struct Burgers1D <: AbstractEquation1D
    max_abs_u0::Float64
end

num_vars(::Burgers1D) = 1

max_wave_speed(eq::Burgers1D) = eq.max_abs_u0

function flux(eq::Burgers1D, Ul::Vector{Float64}, Ur::Vector{Float64})
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

function exact_solution!(utrue::Matrix{Float64}, eq::Burgers1D, xmid::Vector{Float64}, ic::Function, bcs::Dict, x0::Float64, x1::Float64, t::Float64)
    #TODO

end
