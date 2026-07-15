import Gmsh: gmsh


function main()
    gmsh.initialize()

    lc = 1.0

    Lx = 800.0
    Ly = 20.0

    p1 = gmsh.model.geo.addPoint( 0, 0,  0, lc)
    p2 = gmsh.model.geo.addPoint(Lx, 0,  0, lc)
    p3 = gmsh.model.geo.addPoint(Lx, Ly, 0, lc)
    p4 = gmsh.model.geo.addPoint( 0, Ly, 0, lc)

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    c1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])

    s1 = gmsh.model.geo.addPlaneSurface([c1])

    gmsh.model.geo.mesh.setTransfiniteCurve.((l1, l3), 801)
    gmsh.model.geo.mesh.setTransfiniteCurve.((l2, l4), 21)

    gmsh.model.geo.mesh.setTransfiniteSurface(s1)
    gmsh.model.geo.mesh.setRecombine(2, s1)

    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [l1], -1, "Bottom")
    gmsh.model.addPhysicalGroup(1, [l2], -1, "Right")
    gmsh.model.addPhysicalGroup(1, [l3], -1, "Top")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "Left")

    gmsh.model.addPhysicalGroup(2, [s1], -1, "Fluid")

    gmsh.model.mesh.generate(2)

    # perturbing the middle row (±1e-6)

    nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes(2, s1)
    # println(nodeTags)
    # println(coord)

    δy = 1e-6
    odd = 0
    even = 0

    for (i, tag) in enumerate(nodeTags)  
        x = coord[3*i-2]
        y = coord[3*i-1]
        z = coord[3*i]

        if y ≈ 10.0 atol = 1e-2 # target the points that are in the middle line horizontally
            if round(Int, x)%2 == 0  # look up their parity
                gmsh.model.mesh.setNode(tag, [x, y+δy, z], [])
                even +=1
            else
                gmsh.model.mesh.setNode(tag, [x, y-δy, z], [])
                odd +=1
            end
        end
    end 

    println("odd : ", odd)
    println("even : ", even)




    gmsh.write("odd_even_decoupling.msh")

    # gmsh.fltk.run()

    gmsh.finalize()
end

main() 

