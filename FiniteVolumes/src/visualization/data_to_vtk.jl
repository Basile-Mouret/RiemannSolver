using WriteVTK

# Build the 3 × npoints node array (z = 0 for a 2D mesh) used by VTK.
function _vtk_points(mesh::Mesh2D)
    points = zeros(3, length(mesh.points))
    @inbounds for (i, p) in enumerate(mesh.points)
        points[1, i] = p[1]
        points[2, i] = p[2]
    end
    return points
end

# Map node count to a VTK cell type.
_vtk_celltype(n) = n == 3 ? VTKCellTypes.VTK_TRIANGLE :
                   n == 4 ? VTKCellTypes.VTK_QUAD :
                   error("unsupported cell with $n nodes")

_vtk_cells(mesh::Mesh2D) = [MeshCell(_vtk_celltype(length(c)), c) for c in mesh.cells]

# 1D mesh: points laid out along x, y = z = 0; cells are 2-node line segments.
function _vtk_points(mesh::Mesh1D)
    points = zeros(3, length(mesh.points))
    @inbounds for (i, x) in enumerate(mesh.points)
        points[1, i] = x
    end
    return points
end

_vtk_cells(mesh::Mesh1D) =
    [MeshCell(VTKCellTypes.VTK_LINE, (c[1], c[2])) for c in mesh.cells]

# Write every field declared by `output_fields(eq)` as cell data on `vtk`.
# `U` is the ncells × nvars matrix of conserved states.
function _write_output_fields!(vtk, eq::AbstractEquation, U::AbstractMatrix)
    ncells = size(U, 1)
    nvars  = num_vars(eq)
    for f in output_fields(eq)
        if f.kind == :scalar
            data = Vector{Float64}(undef, ncells)
            @inbounds for i in 1:ncells
                data[i] = f.extract(SVector{nvars}(@view U[i, :]))
            end
            vtk[f.name, VTKCellData()] = data
        elseif f.kind == :vector
            d = length(f.extract(SVector{nvars}(@view U[1, :])))
            # WriteVTK wants components along the first dim: d × ncells
            data = Matrix{Float64}(undef, d, ncells)
            @inbounds for i in 1:ncells
                v = f.extract(SVector{nvars}(@view U[i, :]))
                @views data[:, i] .= v
            end
            vtk[f.name, VTKCellData()] = data
        else
            error("unknown output field kind :$(f.kind) for \"$(f.name)\"")
        end
    end
    return nothing
end

"""
    data_to_vtk(mesh, eq, U_hist, dt_hist; outdir, filename)

Write a whole history of states (one matrix per frame) to a ParaView collection.
Fields are taken from `output_fields(eq)`, so this works for any equation.
"""
function data_to_vtk(
    mesh::AbstractMesh,
    eq::AbstractEquation,
    U_hist::Vector{Matrix{Float64}},
    dt_hist::Vector{Float64};
    outdir::AbstractString = "out",
    filename::AbstractString = "simulation",
)
    timesteps = cumsum(dt_hist)
    @assert length(timesteps) == length(U_hist) "got $(length(U_hist)) states but $(length(timesteps)) times"

    points = _vtk_points(mesh)
    cells  = _vtk_cells(mesh)

    mkpath(outdir)

    paraview_collection(joinpath(outdir, filename)) do pvd
        for (i, t) in enumerate(timesteps)
            fname = joinpath(outdir, "$(filename)_$(lpad(i, 4, '0'))")
            vtk_grid(fname, points, cells) do vtk
                _write_output_fields!(vtk, eq, U_hist[i])
                pvd[t] = vtk
            end
        end
    end
end
