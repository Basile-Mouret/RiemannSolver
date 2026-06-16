# ─────────────────────────────────────────────────────────────────────────────
# mesh_generation_example.jl
#
# Conforming 2D mesh with a HOLE and an INCLUSION using Gmsh.jl (OCC kernel).
#
# Geometry
#   Domain    :  rectangle [0, Lx] × [0, Ly]
#   Hole      :  circle at (0.5, 0.5), r=0.15  → void, removed from mesh
#   Inclusion :  circle at (1.5, 0.5), r=0.20  → separate material tag
#
# Physical groups exported
#   dim=2  "Matrix"             – bulk elements
#   dim=2  "Inclusion"          – inclusion elements
#   dim=1  "Left","Right","Top","Bottom"  – outer BC edges
#   dim=1  "HoleBoundary"       – inner wall of the hole
#   dim=1  "InclusionInterface" – shared matrix/inclusion interface
#
# Install
#   julia> import Pkg; Pkg.add("Gmsh")
#   (Gmsh.jl bundles its own Gmsh library — no separate Gmsh install needed)
# ─────────────────────────────────────────────────────────────────────────────

using Gmsh


# ── Geometry parameters ───────────────────────────────────────────────────────

const Lx, Ly = 2.0, 1.0           # domain size

const H_CX, H_CY, H_R = 0.5, 0.5, 0.15    # hole  (removed void)
const I_CX, I_CY, I_R = 1.5, 0.5, 0.20    # inclusion (kept, different tag)


# ── Mesh-size parameters ──────────────────────────────────────────────────────

const LC_FAR    = 0.08   # coarse background element size
const LC_NEAR   = 0.02   # fine element size right at the interfaces
const DIST_RAMP = 0.06   # distance over which the size ramps from fine→coarse


# ─────────────────────────────────────────────────────────────────────────────
# Main meshing function
# ─────────────────────────────────────────────────────────────────────────────

function build_mesh(outfile = "mesh.msh")

    gmsh.initialize()
    gmsh.model.add("domain")
    gmsh.option.setNumber("General.Verbosity", 2)   # 0=silent 5=debug


    # ── 1.  Primitives (OCC kernel) ───────────────────────────────────────────
    #
    # Use the OCC (OpenCASCADE) kernel whenever you need boolean operations.
    # The simpler "geo" kernel does not support fragment / cut / fuse.

    rect = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, Lx, Ly)
    hole = gmsh.model.occ.addDisk(H_CX, H_CY, 0.0, H_R, H_R)
    incl = gmsh.model.occ.addDisk(I_CX, I_CY, 0.0, I_R, I_R)

    # Other useful OCC primitives:
    #   addDisk(cx,cy,0, rx,ry)           → ellipse  (rx ≠ ry)
    #   addRectangle(x,y,0, w,h, r=0)     → rounded rectangle  (r = corner radius)
    #   addBSpline / addSpline + addWire + addSurfaceFilling → arbitrary smooth boundary
    # For a polygon: see the helper at the bottom of this file.


    # ── 2.  Boolean fragment ──────────────────────────────────────────────────
    #
    # fragment(objects, tools) splits all shapes into non-overlapping pieces
    # and creates SHARED boundaries where pieces touch.
    # → The resulting mesh is *conforming* (nodes match at interfaces).
    #
    # Return value:
    #   outMap[1]  → output (dim,tag) pairs that came from `rect`
    #   outMap[2]  → output (dim,tag) pairs that came from `hole` disk
    #   outMap[3]  → output (dim,tag) pairs that came from `incl` disk

    _, outMap = gmsh.model.occ.fragment([(2, rect)], [(2, hole), (2, incl)])

    # !! Always call synchronize() after OCC operations before using the model
    gmsh.model.occ.synchronize()


    # ── 3.  Identify surfaces from the fragment map ───────────────────────────

    hole_surfs = [t for (_, t) in outMap[2]]
    incl_surfs = [t for (_, t) in outMap[3]]
    hole_set   = Set(hole_surfs)
    incl_set   = Set(incl_surfs)
    mat_surfs  = [t for (_, t) in gmsh.model.getEntities(2)
                     if t ∉ hole_set && t ∉ incl_set]

    @info "Surfaces identified" matrix=mat_surfs inclusion=incl_surfs hole=hole_surfs


    # ── 4.  Physical groups – surfaces ────────────────────────────────────────
    #
    # KEY RULE: once any physical group is defined, Gmsh exports ONLY the
    # entities that belong to a group. We intentionally omit hole_surfs
    # so the hole becomes a void in the exported mesh.

    pg = gmsh.model.addPhysicalGroup(2, mat_surfs)
    gmsh.model.setPhysicalName(2, pg, "Matrix")

    pg = gmsh.model.addPhysicalGroup(2, incl_surfs)
    gmsh.model.setPhysicalName(2, pg, "Inclusion")


    # ── 5.  Classify boundary curves ─────────────────────────────────────────
    #
    # getBoundary(dimTags, oriented, combined, recursive)
    #   combined=true → shared curves between surfaces in the list cancel out,
    #                   giving only the outer boundary of their union.

    function bnd_curves(surfs)
        tags = gmsh.model.getBoundary(
            [(2, s) for s in surfs],
            false,   # oriented  (false → return positive tags)
            true,    # combined  (true  → shared curves cancel)
            false    # recursive
        )
        return [abs(t) for (_, t) in tags]
    end

    hole_crvs = bnd_curves(hole_surfs)
    incl_crvs = bnd_curves(incl_surfs)   # inclusion ↔ matrix interface

    # Boundary of (matrix ∪ inclusion): inclusion interface cancels,
    # leaving outer rectangle edges + hole boundary.
    all_outer     = bnd_curves([mat_surfs; incl_surfs])
    hole_crv_set  = Set(hole_crvs)
    rect_crvs     = [t for t in all_outer if t ∉ hole_crv_set]

    # Classify rectangle edges by position (for side-specific BCs)
    left_c = Int[]; right_c = Int[]; top_c = Int[]; bot_c = Int[]
    for t in rect_crvs
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(1, t)
        if abs(xmax - xmin) < 1e-9        # vertical edge
            xmin < 1e-9 ? push!(left_c,  t) : push!(right_c, t)
        else                               # horizontal edge
            ymin < 1e-9 ? push!(bot_c,   t) : push!(top_c,   t)
        end
    end


    # ── 6.  Physical groups – curves ──────────────────────────────────────────

    for (crvs, name) in [
            (left_c,    "Left"),
            (right_c,   "Right"),
            (top_c,     "Top"),
            (bot_c,     "Bottom"),
            (hole_crvs, "HoleBoundary"),
            (incl_crvs, "InclusionInterface"),
        ]
        isempty(crvs) && continue
        pg = gmsh.model.addPhysicalGroup(1, crvs)
        gmsh.model.setPhysicalName(1, pg, name)
    end


    # ── 7.  Mesh-size field ───────────────────────────────────────────────────
    #
    # Distance field  → measures distance from the named curves.
    # Threshold field → maps that distance to an element size.
    # At the interface:        LC_NEAR
    # Far from the interface:  LC_FAR
    # Linear ramp over:        DIST_RAMP

    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "CurvesList",
                                     unique([hole_crvs; incl_crvs]))
    gmsh.model.mesh.field.setNumber(f_dist, "Sampling", 100)

    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField",  f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin",  LC_NEAR)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax",  LC_FAR)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin",  0.0)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax",  DIST_RAMP)

    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

    # Prevent built-in size heuristics from fighting the field
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints",          0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature",       0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary",  0)


    # ── 8.  Meshing algorithm ─────────────────────────────────────────────────
    #
    # 2D algorithms:  1=MeshAdapt  5=Delaunay  6=Frontal  8=Frontal-Delaunay
    # Algorithm 6 (Frontal) produces well-shaped triangles near curved boundaries.
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Laplace2D")   # smooth element quality


    # ── 9.  Export ────────────────────────────────────────────────────────────
    #
    # .msh v4.1 is the right format for modern Julia FEM packages.
    # Use v2.2 only if a legacy package requires it.
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(outfile)
    println("✓ Mesh written → $outfile")

    # Uncomment to inspect the mesh in the Gmsh GUI before saving:
    # gmsh.fltk.run()

    gmsh.finalize()
end


# ─────────────────────────────────────────────────────────────────────────────
# Run
# ─────────────────────────────────────────────────────────────────────────────

build_mesh("mesh.msh")
#   gmsh.write("mesh.vtk")   # ParaView / WriteVTK.jl


# ─────────────────────────────────────────────────────────────────────────────
# Reading the mesh in Julia FEM packages
# ─────────────────────────────────────────────────────────────────────────────
#
# ── Ferrite.jl ──────────────────────────────────────────────────────────────
#   using Ferrite, FerriteMeshParser
#   grid = togrid("mesh.msh")
#   # Access physical groups:
#   addcellset!(grid, "matrix",    "Matrix")
#   addcellset!(grid, "inclusion", "Inclusion")
#   addfacetset!(grid, "left", "Left")
#
# ── GridapGmsh.jl ───────────────────────────────────────────────────────────
#   using Gridap, GridapGmsh
#   model       = GmshDiscreteModel("mesh.msh")
#   Ω_mat       = Triangulation(model, tags=["Matrix"])
#   Ω_inc       = Triangulation(model, tags=["Inclusion"])
#   Γ_left      = BoundaryTriangulation(model, tags=["Left"])
#   Γ_interface = InterfaceTriangulation(model, ["Matrix"], ["Inclusion"])
#
# ── GmshReader.jl (low-level) ────────────────────────────────────────────────
#   using GmshReader
#   msh = GmshReader.load("mesh.msh")
#   # msh.nodes, msh.elements, msh.physicalGroups
#
# ── Other export formats ─────────────────────────────────────────────────────
#   gmsh.write("mesh.vtk")   # ParaView / WriteVTK.jl
#   gmsh.write("mesh.med")   # Code_Aster / Salome
#   gmsh.write("mesh.su2")   # SU2


# ─────────────────────────────────────────────────────────────────────────────
# COMMON VARIATIONS
# ─────────────────────────────────────────────────────────────────────────────

# ── A. Multiple holes / inclusions ───────────────────────────────────────────
#
# function build_multi_mesh()
#     gmsh.initialize(); gmsh.model.add("multi")
#
#     rect = gmsh.model.occ.addRectangle(0, 0, 0, 3.0, 1.0)
#
#     holes = [gmsh.model.occ.addDisk(0.5+k, 0.5, 0, 0.12, 0.12) for k in 0:1]
#     incls = [gmsh.model.occ.addDisk(2.5,   0.5, 0, 0.18, 0.18)]
#
#     all_tools = [(2, t) for t in [holes; incls]]
#     _, outMap = gmsh.model.occ.fragment([(2, rect)], all_tools)
#     gmsh.model.occ.synchronize()
#
#     # outMap[2..1+length(holes)] → hole surfaces
#     # outMap[2+length(holes)..end] → inclusion surfaces
#     hole_surfs = vcat([outMap[1+k] for k in 1:length(holes)]...)  |> tags
#     incl_surfs = vcat([outMap[1+length(holes)+k] for k in 1:length(incls)]...) |> tags
#     # …same pattern as before
# end


# ── B. Polygon inclusion ──────────────────────────────────────────────────────
#
# function make_polygon(xs, ys)
#     pts = [gmsh.model.occ.addPoint(x, y, 0) for (x,y) in zip(xs, ys)]
#     n   = length(pts)
#     lns = [gmsh.model.occ.addLine(pts[i], pts[mod1(i+1, n)]) for i in 1:n]
#     loop = gmsh.model.occ.addCurveLoop(lns)
#     return gmsh.model.occ.addPlaneSurface([loop])
# end
# # Example: hexagonal inclusion
# θ = range(0, 2π, length=7)[1:end-1]
# hex = make_polygon(0.3cos.(θ) .+ 1.5, 0.3sin.(θ) .+ 0.5)


# ── C. Quadrilateral / structured mesh ───────────────────────────────────────
#
# For a quad mesh, call setRecombine after generating:
#   gmsh.model.mesh.generate(2)
#   for (_, t) in gmsh.model.getEntities(2)
#       gmsh.model.mesh.setRecombine(2, t)
#   end
#   gmsh.model.mesh.generate(2)   # re-generate with quads
#
# For a structured transfinite mesh on a rectangle:
#   gmsh.model.mesh.setTransfiniteCurve(tag, nPoints)  # per edge
#   gmsh.model.mesh.setTransfiniteSurface(tag)
#   gmsh.model.mesh.setRecombine(2, tag)


# ── D. Combine multiple size fields ──────────────────────────────────────────
#
# If you have multiple Distance/Threshold pairs, combine with a Min field:
#   f_min = gmsh.model.mesh.field.add("Min")
#   gmsh.model.mesh.field.setNumbers(f_min, "FieldsList", [f_thr1, f_thr2])
#   gmsh.model.mesh.field.setAsBackgroundMesh(f_min)
