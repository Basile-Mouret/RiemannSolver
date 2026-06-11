using LinearAlgebra

struct Advection2D <: AbstractEquation2D
    c::NTuple{2, Float64}
end

Advection2D(; c::NTuple{2, Float64}) = Advection2D(c)

num_vars(::Advection2D) = 2

function flux(eq::Advection2D, UL::Vector{Float64}, UR::Vector{Float64}, normal::NTuple{2, Float64})
    cn = eq.c ⋅ normal
    if cn > 0
        return cn*UL
    else
        return cn*UR
    end
end

function compute_dt(mesh::Mesh2D, eq::Advection2D, values::Matrix{Float64}, CFL::Float64)::Float64
    min_dt = typemax(Float64)
    for cell_idx in 1:length(mesh.cells)
        sum_outgoing = 0.0
        for face_idx in 1:length(mesh.faces)
            if eq.c ⋅ mesh.face_normals[face_idx] > 0
                sum_outgoing += (mesh.face_lengths[face_idx]  * (eq.c ⋅ mesh.face_normals[face_idx]))/mesh.cell_measure[cell_idx]
            end
        end
        dt = CFL / sum_outgoing
        min_dt = min(min_dt, dt)
    end
        
    return min_dt 
end


