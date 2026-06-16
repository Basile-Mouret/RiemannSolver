function explicit_euler_step!(
    new_values::Matrix{Float64},
    cell_values::Matrix{Float64},
    mesh::AbstractMesh,
    eq::AbstractEquation,
    bcs::Dict{String, <: AbstractBC},
    dt::Float64,
    t::Float64
)
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
    for (tag, boundary_faces) in mesh.boundary_tags
        bc = bcs[tag]
        for face_id in boundary_faces
            (CL, CR) = mesh.face_cells[face_id]
            _apply_boundary_face!(new_values, cell_values, mesh, eq, bc, face_id, CL, CR, dt, t)
        end
    end
end

"""
for type stability inside explicit euler step (as the bc is not known at compile time)
"""
function _apply_boundary_face!(
    new_values::Matrix{Float64},
    cell_values::Matrix{Float64},
    mesh::AbstractMesh,
    eq::AbstractEquation,
    bc::AbstractBC,
    face_id::Int,
    CL::Int,
    CR::Int,
    dt::Float64,
    t::Float64
)
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
