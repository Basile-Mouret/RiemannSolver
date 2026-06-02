- TP's
following the lab works, I implemented 1D finite volume solvers for the advection equation and the wave system.
- I created a structure Mesh1D that holds all my meshing information (points, cells, interfaces etc.)
- I created plotting functions, where I can also add the true solution and an animation method

- saved all the data into a U_hist Vector/Matrix => then preprocessing for analysis (Energy, animation, initial conditions etc.)

- I now wanted to refactor the code such that I can reuse the same mesh and solver structures for other types of equations (Burgers and then euler system)

-> DeepSeek V4 Pro with Opencode : cost 17cts in 30 min AI work and 1h to read and test everything

- single solve entry point for all 1D equations
- BCs in a dictionnary
- consistent types (Matrix{Float64}) between scalar and vector
- Initial conditions returns a vector -> uniformity between scalar and vector
- organized into multiple folders and files inside of the FiniteVolumes.jl Module

-> what I would not have done : the abstract types for boundary conditions and equations

No comments in the code



A quel point l'IA est elle utilisée ici?

-> pour débugger? ex: échanger u et p dans mes calculs

9h45 début commenter code

- mettre U_hist_exact en optionel dans l'animation

When we had boundaries I originaly directly computed the flow. Here the LLM replaced it with ghosts cells. For simple Outflow, we set U ghost = U interior that creates no gradient and such the fluid can escape the boundaries. 
