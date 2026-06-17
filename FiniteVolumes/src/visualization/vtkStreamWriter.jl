using WriteVTK

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
    VTKStreamWriter(eq, pvd, points, cells, outdir, filename, Float64(dt_out), 0.0, 0)
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
