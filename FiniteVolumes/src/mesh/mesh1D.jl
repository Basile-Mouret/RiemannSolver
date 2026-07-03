"""
1D Mesh structure

- `points`        : vertex coordinates (only x for 1D)
- `cells`         : vertex indices per cell
- `cell_centers`  : cell barycenters 
- `cell_measure`  : cell length
- `faces`         : vertex index defining a face
- `face_centers`  : face midpoints (equivalent to face vertex for 1D)
- `face_lengths`  : face lengths (for mutli D compat, equals 1 for 1D)
- `face_cells`    : (left, right) cell indices, 0 means exterior (boundary side)
- `face_normals`  : unit outward normal from the left cell at each face
- `boundary_faces`: boundary faces
- `boundary_tags` : correspondance between physical groups and faces
- `cell_perimeters`: sum of face lengths bordering each cell
"""
struct Mesh1D <: AbstractMesh
    points         :: Vector{Float64}
    cells          :: Vector{SVector{2,Int}}
    cell_centers   :: Vector{Float64}
    cell_measure   :: Vector{Float64}
    faces          :: Vector{Int}
    face_centers   :: Vector{Float64}
    face_lengths   :: Vector{Float64}
    face_cells     :: Vector{SVector{2,Int}}
    face_normals   :: Vector{Vector{Float64}}
    boundary_faces :: Vector{Int}
    boundary_tags  :: Dict{String, Vector{Int}}
end

"""
1D regular mesh generation
"""
function generate_1DMesh(x0::Float64, x1::Float64, N::Int, periodic::Bool)::Mesh1D
    points       = collect(range(x0, x1, length = N + 1))
    cells        = [(i, i + 1) for i in 1:N]
    cell_centers = (points[1:end-1] .+ points[2:end]) .* 0.5
    cell_measure = diff(points)

    if periodic
        faces          = collect(1:N)
        face_cells     = SVector{2,Int}[(i, mod(i, N) + 1) for i in 1:N]
        face_normals   = [[1.0] for _ in 1:N]
        boundary_faces = Int[]
        boundary_tags  = Dict{String, Vector{Int}}()
    else
        faces          = collect(1:N+1)
        face_cells     = SVector{2,Int}[(0, 1);[(i, i + 1) for i in 1:(N - 1)];(N, 0)]
        face_normals   = [[1.0] for _ in 1:N+1]
        boundary_faces = [1, N + 1]
        boundary_tags  = Dict("left" => [1], "right" => [N + 1])
    end

    face_centers = points[faces]
    face_lengths = ones(Float64, length(faces))

    return Mesh1D(
                  points,
                  cells,
                  cell_centers,
                  cell_measure,
                  faces,
                  face_centers,
                  face_lengths,
                  face_cells,
                  face_normals,
                  boundary_faces,
                  boundary_tags
                 )
end

"""
loads a 1D mesh from a gmsh file
"""
function load_mesh1D(filename::String)::Mesh1D
    #TODO
end

"""
computes 1D second order quadrature 
"""
function quadrature_1D(x1::Float64, x2::Float64, f::Function)
    h = x2 - x1
    return h / 2 * (f(-h / (2 * sqrt(3)) + (x1 + x2) / 2) + f(h / (2 * sqrt(3)) + (x1 + x2) / 2))
end

"""
computes L2 error between an approximated and a true solution on a 1D mesh
"""
function compute_L2_1D(mesh::Mesh1D, u::Vector{Float64}, u_exact::Function)
    error_val = 0.0
    for (i, (n1, n2)) in enumerate(mesh.cells)
        x1 = mesh.points[n1]
        x2 = mesh.points[n2]
        error_val += quadrature_1D(x1, x2, x -> (u[i] - u_exact(x))^2)
    end
    return sqrt(error_val)
end
