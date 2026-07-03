lc = 0.1;

H = 3.0;
h = 1.0;
L = 3.0;
l = 1.5;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {L,0.0,0.0,lc};
Point(3) = {0.0,H,0.0,lc};
Point(4) = {-L,0.0,0.0,lc};
Point(5) = {-l,0.0,0.0,lc};
Point(6) = {0.0,h,0.0,lc};
Point(7) = {l,0.0,0.0,lc};


Line(1) = {4,5};

Ellipse(2) = {5,1,7,6};
Ellipse(3) = {6,1,7,7};

Line(4) = {7,2};

Ellipse(5) = {2,1,4,3};
Ellipse(6) = {3,1,4,4};

Line Loop(1) = {1, 2, 3, 4, 5, 6};

Plane Surface(1) = {1};


// physical groups

Physical Curve("Bottom") = {1, 4};
Physical Curve("Top") = {5,6};
Physical Curve("Object") = {2,3};

Physical Surface("Domain") = {1};

