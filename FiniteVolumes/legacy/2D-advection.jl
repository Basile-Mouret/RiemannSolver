using CairoMakie, GeometryBasics, LinearAlgebra, Observables
include("mesh2D.jl")

function show_heatmap(mesh::Mesh2D, cellValues::Vector{Float64}, title::String)
    fig = Figure(fontsize=18, size=(600, 600))
    ax = Axis(fig[1, 1], title = title, aspect = DataAspect(), xlabel = "x", ylabel = "y")

    # 1. Build a discrete Polygon for each cell
    # This explicitly maps the exact coordinates to a closed shape, 
    # preventing Makie from scrambling or sharing vertices.
    polys = [
        Polygon([
            Point2f(mesh.points[c[1]]), 
            Point2f(mesh.points[c[2]]), 
            Point2f(mesh.points[c[3]])
        ]) 
        for c in mesh.cells
    ]
    # 2. Plot the polygons
    # Setting strokewidth = 0 removes the triangle borders for a smooth heatmap look
    m = poly!(ax, polys, color = cellValues, colormap = :viridis, strokewidth = 1)

    # Add Colorbar
    Colorbar(fig[1, 2], m, label="Cell Values")
    
    return fig
end


"""
computes 2D transport flux
"""
function flux_transport2D(c::Tuple{Float64,Float64}, U1::Float64, U2::Float64, normal::Tuple{Float64,Float64})
    cn = c ⋅ normal
    if cn > 0
        return cn*U1
    else
        return cn*U2
    end
end


"""
solves the transport equation on a 2D mesh
"""
function solve_transport_2D(mesh::Mesh2D, c::Tuple{Float64,Float64}, u0::Function, n_timesteps::Int, dt::Float64)
    # u_t + cxu_x + cy u_y = 0

    # compute the values for each cell (center)
    cell_values = zeros(Float64, length(mesh.cells))
    for cell_id in 1:length(mesh.cells)
        cell_values[cell_id] = u0(mesh.cells_center[cell_id]...)
    end

    display(show_heatmap(mesh, cell_values, "Initial Condition"))

    color_obs = Observable(cell_values)
    fig = Figure(fontsize=18, size=(600, 600))
    ax = Axis(fig[1, 1], title = "2D Transport", aspect = DataAspect(), xlabel = "x", ylabel = "y")

    polys = [
        Polygon([
            Point2f(mesh.points[u]), 
            Point2f(mesh.points[v]), 
            Point2f(mesh.points[w])
        ]) 
        for (u,v,w) in mesh.cells
    ]

    m = poly!(ax, polys, color = color_obs, colormap = :viridis, strokewidth = 1, colorrange = (-1.0, 1.0))
    Colorbar(fig[1, 2], m, label="Cell Values")

    # 4. The Animation Loop
    # `record` handles the loop and saves the MP4 frame-by-frame
    record(fig, "transport_animation.mp4", 1:n_timesteps; framerate = 30) do step
        newCellValues = copy(cell_values)
        
        for faceIndex in 1:length(mesh.faces)
                cell1 = mesh.face_cells[faceIndex][1]
                cell2 = mesh.face_cells[faceIndex][2]
            if cell1 != 0 && cell2 != 0
                
                u1 = cell_values[cell1]
                u2 = cell_values[cell2]

                # compute the normal vector from cell cell1 to cell cell2
                x1, y1 = mesh.points[mesh.faces[faceIndex][1]]
                x2, y2 = mesh.points[mesh.faces[faceIndex][2]]
                normal = (y2-y1, x1-x2)

                A1 = mesh.cells_area[cell1]
                A2 = mesh.cells_area[cell2]

                F = flux_transport2D(c, u1, u2, normal)
                newCellValues[cell1] -= dt*F/A1
                newCellValues[cell2] += dt*F/A2
            end
        end
        
        cell_values .= newCellValues
        xtrue =  getindex.(mesh.cells_center, 1) .- c[1] .* step.* dt
        ytrue =  getindex.(mesh.cells_center, 2) .- c[2] .* step.* dt
        utrue = u0.(xtrue, ytrue)
        l2_error = sqrt(sum((cell_values .- utrue).^2) / length(cell_values))
        println("Time step $step: L2 error = $l2_error")
        color_obs[] = cell_values 
        ax.title = "Time step $step"
    end
end


function main()
    x0, x1, y0, y1 = 0, 1, 0, 1
    triangulation = generate_mesh2D_rect(x0, x1, y0, y1; seed=67)
    mesh = extract_mesh2D_data(triangulation)
    
    u0(x,y) = sin(2π*x)*cos(2π*y)
    n_timesteps = 10
    dt = 0.01
    c = (1., 2.)

    # TODO boundary conditions
    solve_transport_2D(mesh, c, u0, n_timesteps, dt)
end

main()
