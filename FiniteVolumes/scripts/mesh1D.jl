"""
1D Mesh
"""
struct Mesh1D
    x::Vector{Float64}
    Cells::Vector{Tuple{Int, Int}}
    Faces::Vector{Tuple{Int, Int}}
    BoundaryFaces::Vector{Int}
end

"""
generates a 1D Mesh
"""
function generate_1DMesh(x0::Float64, x1::Float64, N::Int, periodic::Bool)
    x = collect(range(x0, x1, length=N+1))
    Cells = [(i,i+1) for i in 1:N]
    Faces = [(i,i+1) for i in 1:N-1]
    if periodic
        push!(Faces, (N, 1))
        BoundaryFaces = Int[]
    else
        BoundaryFaces = [1, N]
    end
    return Mesh1D(x, Cells, Faces, BoundaryFaces)
end

