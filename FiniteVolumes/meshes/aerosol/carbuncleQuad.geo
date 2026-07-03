R=1;        // Outer boundary radius
r=0.5;        // Cylinder radius
nR = 80;      // Radial discretisation
nC = 80;      // Demi-cirumference discretisation (multiplied by 2)
lc1 = 0.05;   // Caracteristic lenght of the optional Domain extension.

// Domain : 2D cylinder in a 2D circular plan
Point(1) = {0,r,0,1};
l[] = Extrude {0,R,0} {Point{1}; Layers{nR};};
s1[] = Extrude {{0,0,1},{0,0,0},Pi/2}{ Line{l[1]}; Layers{nC}; Recombine; };
s2[] = Extrude {{0,0,1},{0,0,0},Pi/2}{ Line{s1[0]}; Layers{nC}; Recombine; };
//s3[] = Extrude {{0,0,1},{0,0,0},Pi/2}{ Line{s2[0]}; Layers{nC}; Recombine; };
//s4[] = Extrude {{0,0,1},{0,0,0},Pi/2}{ Line{s3[0]}; Layers{nC}; Recombine; };

// NB : for triangular meshes instead of quad remove "Recombine;" from the lines above.


Physical Curve("Outlet") = {l[1], s2[0]};
Physical Curve("Inlet") = {s1[2], s2[2]};
Physical Curve("Cylinder") = {s1[3], s2[3]};

Physical Surface("Domain") = {s1[1], s2[1]};
