#ENV["GKSwstype"] = "100" # no window
using Plots

u0(x) = abs(x - 0.5) < 0.1 ? 1.0 : 0.0



function main()
    N = 100 # number of points
    x0, x1 = 0.0, 1.0 # interval
    dx = (x1 - x0) / N
    c = 1.0 # celerity
    T = 120
    dt =  1/T

    CFL = c*dt/dx
    
    # Grid
    x = range(x0, x1, length=N+1)
    u = u0.(x)
    unew = similar(u)
    

    anim = @animate for step in 1:T
        # Periodic BC's
        @views unew[2:end] .= u[2:end] .- CFL .* (u[2:end] .- u[1:end-1])
        unew[1] = u[1] - CFL * (u[1]-u[N])
        u .= unew
        
        plot(x, u, 
             xlims=(x0, x1), 
             ylims=(-0.5, 2),
             seriestype = :steppost)
    end
    
    
    #mp4(anim, "advection.mp4", fps=60)
    
end

main()
