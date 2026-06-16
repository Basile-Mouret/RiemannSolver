struct Wave2D <: AbstractEquation2D
    kappa::Float64
    rho::Float64
end

Wave2D(; kappa::Float64, rho::Float64) = Wave2D(kappa, rho)

num_vars(::Wave2D) = 3

function flux(eq::Wave2D, uL::AbstractVector{Float64}, uR::AbstractVector{Float64}, normal::NTuple{2, Float64})
    # Godunov's Flux presented in the UPPA course
    pL, ULx, ULy = uL[1], uL[2], uL[3]
    pR, URx, URy = uR[1], uR[2], uR[3]
    c = sqrt(eq.kappa / eq.rho)
    UL_dot_n = ULx * normal[1] + ULy * normal[2]
    UR_dot_n = URx * normal[1] + URy * normal[2]
    Fp = (0.5 / eq.rho) * (UR_dot_n + UL_dot_n) - 0.5 * c * (pR - pL)
    Fu = (0.5 * eq.kappa * (pR + pL) - 0.5 * c * (UR_dot_n - UL_dot_n)) * normal[1]
    Fv = (0.5 * eq.kappa * (pR + pL) - 0.5 * c * (UR_dot_n - UL_dot_n)) * normal[2]
    return SVector(Fp, Fu, Fv)
end

function apply_ghost(bc::Reflecting2D, u_interior::AbstractVector{Float64}, x, ::Float64)
    return SVector(u_interior[1], -u_interior[2], -u_interior[3])
end

function compute_dt(mesh::Mesh2D, eq::Wave2D, values::Matrix{Float64}, CFL::Float64)::Float64
    cell_perimeters = zeros(length(mesh.cells))
    # TODO : should be precomputed in Mesh ...
    for (face_id, (CL, CR)) in enumerate(mesh.face_cells)
        l  = mesh.face_lengths[face_id]
        if CL!=0
            cell_perimeters[CL] += l
        end

        if CR!=0
            cell_perimeters[CR] += l
        end
    end
    return CFL * (2.0 / sqrt(eq.kappa/eq.rho)) * minimum(mesh.cell_measure ./ cell_perimeters)
end

