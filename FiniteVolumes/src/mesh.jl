"""
1D Mesh structure
"""
struct Mesh1D
    x::Vector{Float64}
    Cells::Vector{Tuple{Int, Int}}
    Faces::Vector{Tuple{Int, Int}}
    BoundaryFaces::Vector{Int}
end

"""
1D regular mesh generation
"""
function generate_1DMesh(x0::Float64, x1::Float64, N::Int, periodic::Bool)
    x = collect(range(x0, x1, length = N + 1))
    Cells = [(i, i + 1) for i in 1:N]
    Faces = [(i, i + 1) for i in 1:(N - 1)]
    if periodic
        push!(Faces, (N, 1))
        BoundaryFaces = Int[]
    else
        push!(Faces, (N, 0))
        push!(Faces, (0, 1))
        BoundaryFaces = [1, N]
    end
    return Mesh1D(x, Cells, Faces, BoundaryFaces)
end

"""
computes the center of the cells for a regular 1D mesh
"""
function cell_centers(mesh::Mesh1D)
    return (mesh.x[1:(end - 1)] + mesh.x[2:end]) ./ 2
end

"""
computes the cell width for regular 1D mesh
"""
function cell_width(mesh::Mesh1D)
    N = length(mesh.x) - 1
    return (mesh.x[end] - mesh.x[1]) / N
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
    for (i, (n1, n2)) in enumerate(mesh.Cells)
        x1 = mesh.x[n1]
        x2 = mesh.x[n2]
        error_val += quadrature_1D(x1, x2, x -> (u[i] - u_exact(x))^2)
    end
    return sqrt(error_val)
end
