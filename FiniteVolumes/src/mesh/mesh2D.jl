import Gmsh: gmsh

"""
2D mesh for finite volume methods, supporting triangular and quadrilateral cells.

Fields
- `points`        : vertex coordinates
- `cells`         : vertex indices per cell (3 nodes for triangles, 4 for quads)
- `cell_centers`  : cell barycenters
- `cell_measure`  : cell areas
- `faces`         : vertex index pair (a, b) per face
- `face_centers`  : face midpoints
- `face_cells`    : (left, right) cell indices; 0 means exterior (boundary side)
- `face_normals`  : unit outward normal from the left cell at each face
- `face_lengths`  : face lengths
- `boundary_faces`: indices into `faces` where right cell == 0
- `boundary_tags` : physical group name → face indices (empty if no physical groups)
"""
struct Mesh2D <: AbstractMesh
    points         :: Vector{NTuple{2, Float64}}
    cells          :: Vector{Vector{Int}}
    cell_centers   :: Vector{NTuple{2, Float64}}
    cell_measure   :: Vector{Float64}
    faces          :: Vector{NTuple{2, Int}}
    face_centers   :: Vector{NTuple{2, Float64}}
    face_cells     :: Vector{NTuple{2, Int}}
    face_normals   :: Vector{NTuple{2, Float64}}
    face_lengths   :: Vector{Float64}
    boundary_faces :: Vector{Int}
    boundary_tags  :: Dict{String, Vector{Int}}
end

# Shoelace formula; works for both triangles and convex/concave polygons
function _polygon_area(pts::Vector{NTuple{2, Float64}}) :: Float64
    n = length(pts)
    A = 0.0
    for i in 1:n
        j = mod1(i + 1, n)
        A += pts[i][1] * pts[j][2] - pts[j][1] * pts[i][2]
    end
    return abs(A) / 2
end

"""
    load_mesh2D(filename) -> Mesh2D

Read a Gmsh .msh file (MSH2 or MSH4) and build a `Mesh2D` with full face
connectivity, outward unit normals, and boundary physical-group tags.
"""
function load_mesh2D(filename::String) :: Mesh2D
    gmsh.initialize()
    gmsh.open(filename)

    # Nodes 
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    n_nodes = length(node_tags)
    points = Vector{NTuple{2, Float64}}(undef, n_nodes)
    node_tag_to_idx = Dict{Int, Int}()
    for i in 1:n_nodes
        tag = node_tags[i]
        x = node_coords[3*(i-1) + 1]
        y = node_coords[3*(i-1) + 2]
        points[i] = (x, y)
        node_tag_to_idx[tag] = i
    end

    # 2D elements (triangles = type 2, quads = type 3)
    elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(2)
    cells = Vector{Vector{Int}}()
    for (et, en_tags) in zip(elem_types, elem_node_tags)
        _, _, _, nodes_per_elem, _, _ = gmsh.model.mesh.getElementProperties(et)
        n_elem = length(en_tags) ÷ nodes_per_elem
        for i in 1:n_elem
            verts = [node_tag_to_idx[en_tags[(i-1)*nodes_per_elem + j]] for j in 1:nodes_per_elem]
            push!(cells, verts)
        end
    end
    n_cells = length(cells)

    # Cell geometry
    cell_centers = Vector{NTuple{2, Float64}}(undef, n_cells)
    cell_measure = Vector{Float64}(undef, n_cells)
    for (c, verts) in enumerate(cells)
        pts = [points[v] for v in verts]
        n = length(pts)
        cell_centers[c] = (sum(p[1] for p in pts) / n, sum(p[2] for p in pts) / n)
        cell_measure[c] = _polygon_area(pts)
    end

    # Face connectivity
    # Map canonical edge key (min,max) → list of (cell_id, directed_edge (a,b))
    edge_dict = Dict{NTuple{2,Int}, Vector{Tuple{Int, NTuple{2,Int}}}}()
    for (c, verts) in enumerate(cells)
        n = length(verts)
        for i in 1:n
            a = verts[i]
            b = verts[mod1(i + 1, n)]
            key = minmax(a, b)
            push!(get!(edge_dict, key, []), (c, (a, b)))
        end
    end

    faces          = NTuple{2,Int}[]
    face_centers   = NTuple{2,Float64}[]
    face_cells_vec = NTuple{2,Int}[]
    face_normals   = NTuple{2,Float64}[]
    face_lengths   = Float64[]
    boundary_faces = Int[]
    face_lookup    = Dict{NTuple{2,Int}, Int}()  # canonical key → face index

    for (key, cell_list) in edge_dict
        face_id = length(faces) + 1
        face_lookup[key] = face_id

        a, b = cell_list[1][2]
        push!(faces, (a, b))

        L = cell_list[1][1]
        R = length(cell_list) == 2 ? cell_list[2][1] : 0
        push!(face_cells_vec, (L, R))

        # Unit outward normal from cell L:
        #   edge direction (dx, dy), rotated +90° gives (-dy, dx)
        #   then check sign against (cell_center_L - face_midpoint)
        pa, pb = points[a], points[b]
        dx, dy = pb[1] - pa[1], pb[2] - pa[2]
        len    = sqrt(dx^2 + dy^2)
        nx, ny = -dy / len, dx / len
        mx, my = (pa[1] + pb[1]) / 2, (pa[2] + pb[2]) / 2
        push!(face_centers, (mx, my))
        cl = cell_centers[L]
        # flip if normal points toward L instead of away from it
        if (cl[1] - mx) * nx + (cl[2] - my) * ny > 0
            nx, ny = -nx, -ny
        end
        push!(face_normals, (nx, ny))
        push!(face_lengths, len)

        R == 0 && push!(boundary_faces, face_id)
    end

    # Boundary physical groups → named face index sets
    boundary_tags = Dict{String, Vector{Int}}()
    for (dim, tag) in gmsh.model.getPhysicalGroups(1)
        name = gmsh.model.getPhysicalName(dim, tag)
        isempty(name) && (name = "boundary_$tag")
        group_faces = Int[]
        for ent in gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            _, _, en_tags = gmsh.model.mesh.getElements(1, ent)
            for et_nodes in en_tags
                n_edges = length(et_nodes) ÷ 2
                for j in 1:n_edges
                    a = node_tag_to_idx[et_nodes[2*(j-1) + 1]]
                    b = node_tag_to_idx[et_nodes[2*(j-1) + 2]]
                    fid = get(face_lookup, minmax(a, b), 0)
                    fid > 0 && push!(group_faces, fid)
                end
            end
        end
        boundary_tags[name] = group_faces
    end

    gmsh.finalize()

    return Mesh2D(
        points, cells, cell_centers, cell_measure,
        faces, face_centers, face_cells_vec, face_normals, face_lengths,
        boundary_faces, boundary_tags
    )
end

"""
    face_outward_normal(mesh, face_id, cell_id) -> NTuple{2,Float64}

Return the outward unit normal from `cell_id` across `face_id`.
Flips the stored normal when `cell_id` is the right cell.
"""
function face_outward_normal(mesh::Mesh2D, face_id::Int, cell_id::Int) :: NTuple{2, Float64}
    L, _ = mesh.face_cells[face_id]
    n = mesh.face_normals[face_id]
    return cell_id == L ? n : (-n[1], -n[2])
end
