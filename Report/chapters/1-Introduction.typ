= Introduction

// Present the problematic of this internship : numerical disturbances in compressible fluid simulation, and the carbuncle. get to why this internship is necessary, and what it will bring. then explain how we are going to get there.

When simulating supersonic flows, viscosity is often neglected as this results in simple predictable behaviours. 
These flows can be solved using computationally friendly numerical methods like finite volumes.
However in specific cases these simple methods are plagued by numerical artifacts.
This internship is focused on understanding these numerical methods and the observed instabilities.
A first objective is to implement a first order finite volume method with different Riemann solvers for the Euler equations.
Afterwards, we will study test cases that triggers grid-aligned numerical instabilities and add them to the AeroSol software @AeroSol.
Finally a curl and divergence preserving method develloped by the team at Cagire will be tested for these instabilities.

This internship was done at the Inria team Cagire in Pau.
The team is working on various high performance flow simulationswich include complex physical properties like being multiphasic, turbulent and multiscale.
They also work on efficient high order numerical methods for complex geometries.

