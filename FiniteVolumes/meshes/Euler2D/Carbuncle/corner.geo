lc = 1.0;
L = 7.0;
lx = 2.0;
ly = 3.0;

Point(1) = {lx,0.0,0.0,lc};
Point(2) = {L,0.0,0.0,lc};
Point(3) = {L,ly,0.0,lc};
Point(4) = {L,L,0.0,lc};
Point(5) = {lx,L,0.0,lc};
Point(6) = {0.0,L,0.0,lc};
Point(7) = {0.0,ly,0.0,lc};
Point(8) = {lx,ly,0.0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};

Line(9) = {3,8};
Line(10) = {5,8};

Curve Loop(1) = {1, 2, 9, 8};
Curve Loop(2) = {-9, 3, 4, 10};
Curve Loop(3) = {7, -10, 5, 6};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

Transfinite Curve{1, 9, 4} = 50;
Transfinite Curve{2, 8} = 30; 
Transfinite Curve{3, 10, 6} = 40; 
Transfinite Curve{5, 7} = 20; 

Transfinite Surface {1, 2, 3}; 
Recombine Surface {1, 2, 3}; 

Physical Curve("Bottom") = {1};
Physical Curve("Top") = {4, 5};
Physical Curve("Right") = {2, 3};
Physical Curve("Left") = {6};
Physical Curve("Object") = {7, 8};

Physical Surface("Fluid") = {1, 2, 3};
