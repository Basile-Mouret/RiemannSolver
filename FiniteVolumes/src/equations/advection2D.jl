using LinearAlgebra

struct Advection2D <: AbstractEquation2D
    c::NTuple{2, Float64}
end

Advection2D(; c::NTuple{2, Float64}) = Advection2D(c)

num_vars(::Advection2D) = 1

function flux(eq::Advection2D, UL::Vector{Float64}, UR::Vector{Float64}, normal::NTuple{2, Float64})
    cn = eq.c ⋅ normal
    if cn > 0
        return cn*UL
    else
        return cn*UR
    end
end

function compute_dt(mesh::Mesh2D, eq::Advection2D, values::Matrix{Float64}, CFL::Float64)::Float64
    sum_outgoing = zeros(length(mesh.cells))
    for (face_id, (CL, CR)) in enumerate(mesh.face_cells)
        cn = eq.c ⋅ mesh.face_normals[face_id]
        l  = mesh.face_lengths[face_id]
        if cn > 0 && CL != 0
            sum_outgoing[CL] += cn * l
        elseif cn < 0 && CR != 0
            sum_outgoing[CR] += (-cn) * l
        end
    end
    return minimum(CFL .* mesh.cell_measure ./ sum_outgoing)
end


