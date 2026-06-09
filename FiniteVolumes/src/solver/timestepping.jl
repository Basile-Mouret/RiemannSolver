function explicit_euler_step!(
    new_values::Matrix{Float64},
    values::Matrix{Float64},
    mesh::Mesh1D,
    eq::AbstractEquation1D,
    bcs,
    dt::Float64,
    t::Float64
)
    nvars = num_vars(eq)

    # get the values for the left and right cells
    # if it is a border, use the ghost cells
    for (CL, CR) in mesh.faces_cells
        if CL != 0 && CR != 0
            uL = values[CL, :]
            uR = values[CR, :]
        elseif CL == 0
            uL = apply_ghost(bcs[:left], values[CR, :], t)
            uR = values[CR, :]
        elseif CR == 0
            uL = values[CL, :]
            uR = apply_ghost(bcs[:right], values[CL, :], t)
        end

        # compute the flux for the given equation
        # we should add normal for flux calculations (1 for 1D and face.normal for 2D/3D)
        F = flux(eq, uL, uR)

        # apply fluxes to non ghost cells
        if CL != 0
            for v in 1:nvars
                new_values[CL, v] -= (dt / mesh.cells_length[CL]) * F[v]
            end
        end
        if CR != 0
            for v in 1:nvars
                new_values[CR, v] += (dt / mesh.cells_length[CR]) * F[v]
            end
        end
    end
end
