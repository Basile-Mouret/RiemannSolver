lc = 0.5;

Lx = 10.0;
Ly = 8.0;
lx = 2.0;
ly = (Lx-lx)*Tan(Pi/6.0);

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {lx,0.0,0.0,lc};
Point(3) = {Lx,ly,0.0,lc};
Point(4) = {Lx,Ly,0.0,lc};
Point(5) = {lx,Ly,0.0,lc};
Point(6) = {0.0,Ly,0.0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {2,5};

Curve Loop(1) = {1, 7, 5, 6};
Curve Loop(2) = {2, 3, 4, -7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Curve{3, 6, 7} = 40;
Transfinite Curve{1, 5} = 10; 
Transfinite Curve{2, 4} = 40; 

Transfinite Surface {1, 2}; 
Recombine Surface {1, 2}; 

Physical Curve("Bottom") = {1};
Physical Curve("Top") = {4, 5};
Physical Curve("Right") = {3};
Physical Curve("Left") = {6};
Physical Curve("Object") = {2};

Physical Surface("Fluid") = {1, 2};
