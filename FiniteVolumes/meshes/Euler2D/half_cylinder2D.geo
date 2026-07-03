lc = 0.01;

R = 3.0;
r = 1.0;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {R,0.0,0.0,lc};
Point(3) = {0.0,R,0.0,lc};
Point(4) = {-R,0.0,0.0,lc};
Point(5) = {-r,0.0,0.0,lc};
Point(6) = {0.0,r,0.0,lc};
Point(7) = {r,0.0,0.0,lc};


Line(1) = {4,5};

Circle(2) = {5,1,6};
Circle(3) = {6,1,7};

Line(4) = {7,2};

Circle(5) = {2,1,3};
Circle(6) = {3,1,4};

Line Loop(1) = {1, 2, 3, 4, 5, 6};

Plane Surface(1) = {1};


// physical groups

Physical Curve("Outlet") = {1, 4};
Physical Curve("Inlet") = {5,6};
Physical Curve("Cylinder") = {2,3};

Physical Surface("Domain") = {1};

