#ENV["GKSwstype"] = "100" # no window
using Plots
gr()

u0(x) = abs(x - 0.5) < 0.1 ? 1.0 : 0.0


function main()
    N = 100 # number of points
    x0, x1 = 0.0, 1.0 # interval
    dx = (x1 - x0) / N 
    c = 1.0 # celerity
    T = 120
    dt =  1/T

    CFL = c*dt/dx
    println("CFL = ",CFL)



    
    # Grid 
    x = range(x0, x1, length=N+1)
    cellFaces = hcat(1:N-1, 2:N)'
    cellFaces[2,N-1] = 1
    faceCells = hcat(1:N, 2:N+1)'
    faceCells[2,N] = 1
    borderFaces = []
    #println(cellFaces)
    #println(faceCells)

    xmid = x[1:end-1].+dx/2
    u = u0.(xmid)
    display(plot(xmid, u, seriestype = :steppost))

    unew = similar(u)
    anim = @animate for step in 1:T
        for face in 1:N
            lCell, rCell = faceCells[:,face]
            # as c>0, the initial values are transported to the right
            # we update the right values using previous information from the left
            unew[rCell] = u[rCell] - CFL*(u[rCell]-u[lCell])
        end

        u .= unew
        utrue = u0.(mod.(xmid.-c*step*dt,1.0))
        plot(xmid, utrue,
             xlims=(x0, x1), 
             ylims=(-0.5, 2),
             seriestype = :steppost)

        
        plot!(xmid, u, 
             xlims=(x0, x1), 
             ylims=(-0.5, 2),
             seriestype = :steppost)
    end
    
    #mp4(anim, "advection.mp4", fps=60)
end


main()
