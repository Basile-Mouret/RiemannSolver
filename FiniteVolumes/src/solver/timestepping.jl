function explicit_euler_step!(
    new_values::Matrix{Float64},
    cell_values::Matrix{Float64},
    mesh::M,
    eq::E,
    bcs::Dict{String, BC},
    dt::Float64,
    t::Float64,
) where {M<:AbstractMesh, E<:AbstractEquation, BC<:AbstractBC}
    nvars = num_vars(eq)

    # get the values for the left and right cells
    # if it is a border, use the ghost cells
    # loop over interior cells
    for (face_id, (CL, CR)) in enumerate(mesh.face_cells)
        if CL != 0 && CR != 0
            uL = SVector{nvars}(@view cell_values[CL, :])
            uR = SVector{nvars}(@view cell_values[CR, :])
            # compute the flux for the given equation
            F = flux(eq, uL, uR, mesh.face_normals[face_id])
            coefL = dt * mesh.face_lengths[face_id] / mesh.cell_measure[CL]
            coefR = dt * mesh.face_lengths[face_id] / mesh.cell_measure[CR]
            for v in 1:nvars
                new_values[CL, v] -= coefL * F[v]
                new_values[CR, v] += coefR * F[v]
            end
        end
    end

    # loop over borders
    # dispatch on the (possibly abstract) bc once per tag, then run the
    # per-face loop inside a function barrier where bc is concrete
    for (tag, boundary_faces) in mesh.boundary_tags
        _apply_boundary_faces!(new_values, cell_values, mesh, eq, bcs[tag], boundary_faces, dt, t)
    end
end

"""
function barrier over the faces of a single boundary tag, specialized on the
concrete bc type so the inner per-face work is allocation-free
"""
function _apply_boundary_faces!(
    new_values::Matrix{Float64},
    cell_values::Matrix{Float64},
    mesh::M,
    eq::E,
    bc::BC,
    boundary_faces::Vector{Int},
    dt::Float64,
    t::Float64,
) where {M<:AbstractMesh, E<:AbstractEquation, BC<:AbstractBC}
    for face_id in boundary_faces
        (CL, CR) = mesh.face_cells[face_id]
        nvars = num_vars(eq)
        if CL == 0
            uR = SVector{nvars}(@view cell_values[CR, :])
            uL = apply_ghost(bc, uR, mesh.face_centers[face_id], t)
            F = flux(eq, uL, uR, mesh.face_normals[face_id])
            coefR = dt * mesh.face_lengths[face_id] / mesh.cell_measure[CR]
            for v in 1:nvars
                new_values[CR, v] += coefR * F[v]
            end
        else
            uL = SVector{nvars}(@view cell_values[CL, :])
            uR = apply_ghost(bc, uL, mesh.face_centers[face_id], t)
            F = flux(eq, uL, uR, mesh.face_normals[face_id])
            coefL = dt * mesh.face_lengths[face_id] / mesh.cell_measure[CL]
            for v in 1:nvars
                new_values[CL, v] -= coefL * F[v]
            end
        end
    end
end

