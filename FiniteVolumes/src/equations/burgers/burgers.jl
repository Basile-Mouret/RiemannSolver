struct Burgers1D <: AbstractEquation1D end

num_vars(::Burgers1D) = 1

function flux(::Burgers1D, Ul::AbstractVector{Float64}, Ur::AbstractVector{Float64})
    ul, ur = Ul[1], Ur[1]
    us = 0.0
    if ul > ur # shock
        S = 0.5*(ul+ur) # shock speed
        if S> 0 # right shock
            return us = ul
        else # left shock
            return us = ur
        end
    else # rarefaction
        if ul >= 0 # right supersonic rarefaction
            return us = ul
        elseif ur <= 0 # left supersonic rarefaction
            return us = ur
        else # transonic rarefaction
            return us = 0.0 
        end
    end
    return 0.5*us*us
end

function compute_dt(mesh::Mesh1D, eq::Burgers1D, values::Matrix{Float64}, CFL::Float64)::Float64
    return minimum(CFL .* mesh.cell_measure ./ abs.(values))
end

output_fields(::Burgers1D) = [OutputField("u", :scalar, U -> U[1])]
