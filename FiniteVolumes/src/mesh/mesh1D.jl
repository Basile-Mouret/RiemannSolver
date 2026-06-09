"""
1D Mesh structure
"""
struct Mesh1D
    points         :: Vector{Float64}
    cells          :: Vector{NTuple{2, Int}}
    cells_center   :: Vector{Float64}
    cells_length   :: Vector{Float64}
    faces          :: Vector{Int}
    faces_cells    :: Vector{NTuple{2, Int}}
    faces_normals  :: Vector{Float64}
    boundary_faces :: Vector{Int}
    boundary_tags  :: Dict{String, Vector{Int}}
end

"""
1D regular mesh generation
"""
function generate_1DMesh(x0::Float64, x1::Float64, N::Int, periodic::Bool)
    points       = collect(range(x0, x1, length = N + 1))
    cells        = [(i, i + 1) for i in 1:N]
    cells_center = (points[1:end-1] .+ points[2:end]) .* 0.5
    cells_length = diff(points)

    if periodic
        faces         = collect(1:N)
        faces_cells   = NTuple{2,Int}[(i, mod(i, N) + 1) for i in 1:N]
        faces_normals = ones(Float64, N)
        boundary_faces = Int[]
        boundary_tags  = Dict{String, Vector{Int}}()
    else
        faces       = collect(1:N+1)
        faces_cells = NTuple{2,Int}[(0, 1);[(i, i + 1) for i in 1:(N - 1)];(N, 0)]
        faces_normals = ones(Float64, N + 1)
        faces_normals[1] = -1.0
        boundary_faces = [1, N + 1]       # fixed: was [1, N]
        boundary_tags  = Dict("left" => [1], "right" => [N + 1])
    end

    return Mesh1D(
        points, cells, cells_center, cells_length,
        faces, faces_cells, faces_normals,
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
