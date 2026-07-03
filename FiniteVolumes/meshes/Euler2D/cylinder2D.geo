lc = 0.01;

// rectangle points
xmin = 0.0;
xmax = 5.0;
ymin = -2.0;
ymax = 2.0;

Point(1) = {xmin,ymin,0.0,lc};
Point(2) = {xmax,ymin,0.0,5.0*lc};
Point(3) = {xmax,ymax,0.0,5.0*lc};
Point(4) = {xmin,ymax,0.0,lc};

xc = 2.0;
yc = 0.0;
r = 0.4;

// circle points
Point(5) = {xc, yc, 0.0, lc}; // center
Point(6) = {xc-r, yc, 0.0, lc};
Point(7) = {xc, yc-r, 0.0, lc};
Point(8) = {xc+r, yc, 0.0, 2.0*lc};
Point(9) = {xc, yc+r, 0.0, lc};


// rectangle lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// circle line
Circle(5) = {6,5,7};
Circle(6) = {7,5,8};
Circle(7) = {8,5,9};
Circle(8) = {9,5,6};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1, 2};


// physical groups

Physical Curve("Bottom") = {1};
Physical Curve("Right") = {2};
Physical Curve("Top") = {3};
Physical Curve("Left") = {4};

Physical Curve("Cylinder") = {5, 6, 7, 8};

Physical Surface("Domain") = {1};

