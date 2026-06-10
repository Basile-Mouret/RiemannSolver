"""
1D Mesh structure
"""
struct Mesh1D <: AbstractMesh
    points         :: Vector{Float64}
    cells          :: Vector{NTuple{2, Int}}
    cell_centers   :: Vector{Float64}
    cell_measure   :: Vector{Float64}
    faces          :: Vector{Int}
    face_cells    :: Vector{NTuple{2, Int}}
    face_normals  :: Vector{Vector{Float64}}
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
        faces         = collect(1:N)
        face_cells    = NTuple{2,Int}[(i, mod(i, N) + 1) for i in 1:N]
        face_normals  = [[1.0] for _ in 1:N+1]
        boundary_faces = Int[]
        boundary_tags  = Dict{String, Vector{Int}}()
    else
        faces        = collect(1:N+1)
        face_cells   = NTuple{2,Int}[(0, 1);[(i, i + 1) for i in 1:(N - 1)];(N, 0)]
        face_normals  = [[1.0] for _ in 1:N+1]
        boundary_faces = [1, N + 1]
        boundary_tags  = Dict("left" => [1], "right" => [N + 1])
    end

    return Mesh1D(
        points, cells, cell_centers, cell_measure,
        faces, face_cells, face_normals,
        boundary_faces, boundary_tags
    )
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
