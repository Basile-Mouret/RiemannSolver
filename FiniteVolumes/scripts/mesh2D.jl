using DelaunayTriangulation, StableRNGs
# TODO: change for GMSH : load a file, convert to custom type (! boundary conditions), generate mesh, save mesh, 
# - remove triangulation in mesh2D
# - use Svectors instead of tuples
"""
2D Mesh
"""
struct Mesh2D
    triangulation::DelaunayTriangulation.Triangulation      # triangulation
    points::Vector{Tuple{Float64, Float64}}                 # vertexes
    cells::Vector{Tuple{Int, Int, Int}}                     # vertex Ids forming the cell
    cells_center::Vector{Tuple{Float64, Float64}}           # cell centers
    cells_area::Vector{Float64}                             # cell Area
    faces::Vector{Tuple{Int, Int}}                          # vertex Ids forming the edge
    face_cells::Vector{Tuple{Int, Int}}                     # left and right cell Ids
    boundary_faces::Vector{Int}                             # indices of faces on the domain boundary
end

"""
generates a 2D Mesh using DelaunayTriangulation.jl
"""
function generate_mesh2D_rect(x0, x1, y0, y1; seed=67)
    boundary_points = [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
    boundary_nodes, points = convert_boundary_points_to_indices(boundary_points)
    tri = triangulate(points; boundary_nodes)
    rng = StableRNG(seed)
    refine!(tri; max_area = 0.001 * get_area(tri), rng = rng) 

    # Return the generated triangulation object
    return tri
end

"""
Extracts data from a triangulation into a Mesh2D struct
"""
function extract_mesh2D_data(tri::DelaunayTriangulation.Triangulation)
    # points
    n_points = DelaunayTriangulation.num_points(tri)
    points = [Tuple(get_point(tri, i)) for i in 1:n_points]
    
    # cells
    cells = Tuple{Int, Int, Int}[]
    cells_centers = Tuple{Float64,Float64}[]
    cells_area = Float64[]
    
    # Dictionary to map a sorted tuple of vertices to a specific Cell ID (1 to N)
    triangle_to_cell_id = Dict{Tuple{Int, Int, Int}, Int}()

    for (cell_id, T) in enumerate(each_solid_triangle(tri))
        u, v, w = T[1], T[2], T[3]
        push!(cells, (u, v, w))
        x1, y1 = points[u]
        x2, y2 = points[v]
        x3, y3 = points[w]

        # centers
        cx = (x1+x2+x3)/3.0
        cy = (y1+y2+y3)/3.0 
        push!(cells_centers,(cx,cy))

        # area
        Area = 0.5 * abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))
        push!(cells_area, Area)
        
        # store sorted vertices
        sorted_vertices = Tuple(sort([u, v, w]))
        triangle_to_cell_id[sorted_vertices] = cell_id
    end
    
    # edge informations
    faces = Tuple{Int, Int}[]
    face_cells = Tuple{Int, Int}[]
    boundary_faces = Int[]
    
    for (face_id, e) in enumerate(each_solid_edge(tri))
        u, v = e[1], e[2]
        push!(faces, (u, v))
        
        # Find Adjacent Vertices
        v_left = get_adjacent(tri, u, v)
        v_right = get_adjacent(tri, v, u)
        
        # Determine Cell IDs (0 if it is a boundary)
        left_cell_id = 0
        if !DelaunayTriangulation.is_ghost_vertex(v_left)
            left_cell_id = triangle_to_cell_id[Tuple(sort([u, v, v_left]))]
        end
        
        right_cell_id = 0
        if !DelaunayTriangulation.is_ghost_vertex(v_right)
            right_cell_id = triangle_to_cell_id[Tuple(sort([v, u, v_right]))]
        end
        
        push!(face_cells, (left_cell_id, right_cell_id))
        
        # tag boundary faces
        if left_cell_id == 0 || right_cell_id == 0
            push!(boundary_faces, face_id)
        end
    end
    
    return Mesh2D(tri, points, cells, cells_centers, cells_area, faces, face_cells, boundary_faces)
end

