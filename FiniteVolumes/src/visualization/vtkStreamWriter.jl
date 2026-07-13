using WriteVTK

"""
build the 3xn node array, sets unused dimensions to 0 (for 1D and 2D) as vtk points are always defined in a 3D space
"""
function _vtk_points(mesh::Mesh1D)
    points = zeros(3, length(mesh.points))
    @inbounds for (i, x) in enumerate(mesh.points)
        points[1, i] = x
    end
    return points
end

function _vtk_points(mesh::Mesh2D)
    points = zeros(3, length(mesh.points))
    @inbounds for (i, p) in enumerate(mesh.points)
        points[1, i] = p[1]
        points[2, i] = p[2]
    end
    return points
end

# 1D mesh: points laid out along x, y = z = 0; cells are 2-node line segments.
_vtk_cells(mesh::Mesh1D) = [MeshCell(VTKCellTypes.VTK_LINE, (c[1], c[2])) for c in mesh.cells]
# finds the 2D cell type corresponding to the number of nodes (only triangles and quad for now)
_vtk_2Dcelltype(n) = n == 3 ? VTKCellTypes.VTK_TRIANGLE :
                   n == 4 ? VTKCellTypes.VTK_QUAD :
                   error("unsupported cell with $n nodes")

_vtk_cells(mesh::Mesh2D) = [MeshCell(_vtk_2Dcelltype(length(c)), c) for c in mesh.cells]

function _write_output_fields!(vtk, eq::AbstractEquation, U::AbstractVector{SVector{N, T}}) where {T<:Real, N}
    for f in output_fields(eq)
        if f.kind == :scalar
            data = f.extract.(U)
            vtk[f.name, VTKCellData()] = data
        elseif f.kind == :vector
            vals = f.extract.(U)                       # Vector{SVector{d,T}}
            vtk[f.name, VTKCellData()] = reinterpret(reshape, T, vals)   # d × ncells view
        else
            error("unknown output field kind :$(f.kind) for \"$(f.name)\"")
        end
    end
    return nothing
end

mutable struct VTKStreamWriter{T<:Real, E<:AbstractEquation}
    eq     :: E
    pvd    :: WriteVTK.CollectionFile
    points :: Matrix{T}            # 3 × npoints, built once
    cells  :: Vector{MeshCell}
    dir    :: String
    name   :: String
    dt_out :: T
    next_t :: T
    frame  :: Int
end

function VTKStreamWriter(mesh::AbstractMesh, eq::AbstractEquation, dt_out::T;
                         outdir::AbstractString = "out",
                         filename::AbstractString = "simulation") where {T<:Real}
    mkpath(outdir)
    points = _vtk_points(mesh)
    cells  = _vtk_cells(mesh)
    pvd    = paraview_collection(joinpath(outdir, filename))
    return VTKStreamWriter{T, typeof(eq)}(eq, pvd, points, cells, String(outdir), String(filename), dt_out, zero(T), 0)
end

# Write the current state as a frame at time t
function write_frame!(w::VTKStreamWriter, U, t::Float64)
    U_host = Array(U)
    w.frame += 1
    fname = joinpath(w.dir, "$(w.name)_$(lpad(w.frame, 5, '0'))")
    vtk_grid(fname, w.points, w.cells; compress=false) do vtk
        _write_output_fields!(vtk, w.eq, U_host)
        w.pvd[t] = vtk            # register frame at time t; vtu saved at block end
    end
    return nothing
end

# Call every step; writes only when simulation time crosses the next interval
function maybe_write!(w::VTKStreamWriter, U, t::Float64)    if t >= w.next_t
        write_frame!(w, U, t)
        # advance, skipping any intervals jumped over by a large dt
        while w.next_t <= t
            w.next_t += w.dt_out
        end
    end
    return nothing
end

close_writer!(w::VTKStreamWriter) = vtk_save(w.pvd)
