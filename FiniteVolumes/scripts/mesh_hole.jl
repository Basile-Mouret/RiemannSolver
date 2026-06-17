using Gmsh

const Lx, Ly          = 2.0, 1.0
const H_CX, H_CY, H_R = 0.5, 0.5, 0.15

const LC_FAR    = 0.08
const LC_NEAR   = 0.02
const DIST_RAMP = 0.06

function build_mesh(outfile = "mesh.msh")

    gmsh.initialize()
    gmsh.model.add("domain")
    gmsh.option.setNumber("General.Verbosity", 2)

    rect = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, Lx, Ly)
    hole = gmsh.model.occ.addDisk(H_CX, H_CY, 0.0, H_R, H_R)
    gmsh.model.occ.cut([(2, rect)], [(2, hole)])
    gmsh.model.occ.synchronize()

    mat_surfs = [t for (_, t) in gmsh.model.getEntities(2)]
    pg = gmsh.model.addPhysicalGroup(2, mat_surfs)
    gmsh.model.setPhysicalName(2, pg, "Domain")

    all_crvs = let tags = gmsh.model.getBoundary([(2, s) for s in mat_surfs], false, true, false)
        [abs(t) for (_, t) in tags]
    end

    # Axis-aligned curves → rectangle sides; curved curves → hole boundary
    left_c = Int[]; right_c = Int[]; top_c = Int[]; bot_c = Int[]; hole_crvs = Int[]
    for t in all_crvs
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(1, t)
        if abs(xmax - xmin) < 1e-6
            xmin < 1e-6 ? push!(left_c, t) : push!(right_c, t)
        elseif abs(ymax - ymin) < 1e-6
            ymin < 1e-6 ? push!(bot_c, t) : push!(top_c, t)
        else
            push!(hole_crvs, t)
        end
    end

    for (crvs, name) in [(left_c, "Left"), (right_c, "Right"),
                         (top_c, "Top"), (bot_c, "Bottom"),
                         (hole_crvs, "HoleBoundary")]
        isempty(crvs) && continue
        pg = gmsh.model.addPhysicalGroup(1, crvs)
        gmsh.model.setPhysicalName(1, pg, name)
    end

    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "CurvesList", hole_crvs)
    gmsh.model.mesh.field.setNumber(f_dist, "Sampling", 100)

    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField",  f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin",  LC_NEAR)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax",  LC_FAR)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin",  0.0)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax",  DIST_RAMP)
    
    gmsh.option.setNumber("Mesh.MeshSizeFactor", 1.0)

    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints",         0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature",      0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Laplace2D")

    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(outfile)
    println("✓ $outfile")

    # gmsh.fltk.run()
    gmsh.finalize()
end

build_mesh("meshes/mesh_hole.msh")
