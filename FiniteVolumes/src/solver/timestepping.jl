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
            uL = cell_values[CL, :]
            uR = cell_values[CR, :]
            # compute the flux for the given equation
            F = flux(eq, uL, uR, mesh.face_normals[face_id])
            for v in 1:nvars
                new_values[CL, v] -= (dt / mesh.cell_measure[CL]) * F[v]
                new_values[CR, v] += (dt / mesh.cell_measure[CR]) * F[v]
            end
        end
    end

    # loop over borders 
    for (tag, boundary_faces) in mesh.boundary_tags
        for face_id in boundary_faces
            (CL, CR) = mesh.face_cells[face_id]
            if CL==0
                uR = cell_values[CR, :]
                uL = apply_ghost(bcs[tag], uR, mesh.face_centers[face_id], t)
                F = flux(eq, uL, uR, mesh.face_normals[face_id])
                for v in 1:nvars
                    new_values[CR, v] += (dt / mesh.cell_measure[CR]) * F[v]
                end
            else
                uL = cell_values[CL, :]
                uR = apply_ghost(bcs[tag], uL, mesh.face_centers[face_id], t)
                F = flux(eq, uL, uR, mesh.face_normals[face_id])
                for v in 1:nvars
                    new_values[CL, v] -= (dt / mesh.cell_measure[CL]) * F[v]
                end
            end
        end
    end
end
