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


function _fill_scalar!(data::Vector{Float64}, extract::F, U::AbstractMatrix, ::Val{N}) where {F, N}
    @inbounds for i in eachindex(data)
        data[i] = extract(SVector{N}(@view U[i, :]))
    end
    return nothing
end

function _fill_vector!(data::Matrix{Float64}, extract::F, U::AbstractMatrix, ::Val{N}) where {F, N}
    @inbounds for i in axes(data, 2)
        v = extract(SVector{N}(@view U[i, :]))
        @views data[:, i] .= v
    end
    return nothing
end

function _write_output_fields!(vtk, eq::AbstractEquation, U::AbstractMatrix)
    ncells = size(U, 1)
    nv     = num_vars(eq)
    nvars  = Val(nv)
    for f in output_fields(eq)
        if f.kind == :scalar
            data = Vector{Float64}(undef, ncells)
            _fill_scalar!(data, f.extract, U, nvars)
            vtk[f.name, VTKCellData()] = data
        elseif f.kind == :vector
            d = length(f.extract(SVector{nv}(@view U[1, :])))   # runtime probe, once per field
            data = Matrix{Float64}(undef, d, ncells)
            _fill_vector!(data, f.extract, U, nvars)
            vtk[f.name, VTKCellData()] = data
        else
            error("unknown output field kind :$(f.kind) for \"$(f.name)\"")
        end
    end
    return nothing
end

mutable struct VTKStreamWriter{E<:AbstractEquation}
    eq       :: E                # equation: provides the output fields
    pvd                          # the open collection handle
    points   :: Matrix{Float64}  # 3 × npoints, built once
    cells    :: Vector{<:MeshCell}
    dir      :: String
    name     :: String
    dt_out   :: Float64          # physical interval between frames
    next_t   :: Float64          # next time we should write
    frame    :: Int              # frame counter (for filenames)
end

function VTKStreamWriter(mesh::AbstractMesh, eq::AbstractEquation;
                         outdir::AbstractString = "out",
                         filename::AbstractString = "simulation",
                         dt_out::Real)
    mkpath(outdir)

    points = _vtk_points(mesh)
    cells  = _vtk_cells(mesh)

    pvd = paraview_collection(joinpath(outdir, filename))  # stays open
    VTKStreamWriter(eq, pvd, points, cells, String(outdir), String(filename), Float64(dt_out), 0.0, 0)
end

# Write the current state as a frame at time t
function write_frame!(w::VTKStreamWriter, U::AbstractMatrix, t::Real)
    w.frame += 1
    fname = joinpath(w.dir, "$(w.name)_$(lpad(w.frame, 5, '0'))")
    vtk_grid(fname, w.points, w.cells; compress=false) do vtk
        _write_output_fields!(vtk, w.eq, U)
        w.pvd[t] = vtk            # register frame at time t; vtu saved at block end
    end
    return nothing
end

# Call every step; writes only when simulation time crosses the next interval
function maybe_write!(w::VTKStreamWriter, U::AbstractMatrix, t::Real)
    if t >= w.next_t
        write_frame!(w, U, t)
        # advance, skipping any intervals jumped over by a large dt
        while w.next_t <= t
            w.next_t += w.dt_out
        end
    end
    return nothing
end

close_writer!(w::VTKStreamWriter) = vtk_save(w.pvd)
