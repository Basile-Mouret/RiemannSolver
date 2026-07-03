lc = 0.1;


// rectangle points
xmin = 0.0;
xmax = 5.0;
ymin = -2.0;
ymax = 2.0;

Point(1) = {xmin,ymin,0.0,lc};
Point(2) = {xmax,ymin,0.0,lc};
Point(3) = {xmax,ymax,0.0,lc};
Point(4) = {xmin,ymax,0.0,lc};

t = 0.12;
N = 50;
x0 = 2;

topPoints[] = {};
botPoints[] = {4+1};

For i In {0:N}
    x = 1 - Cos((Pi/2) * (i/N));
    y = 5.0*t*(0.2969*Sqrt(x) - 0.1260*x - 0.3516*x*x + 0.2843*x*x*x - 0.1015*x*x*x*x);

    // top
    Point(4 + i + 1) = {x0+x, y, 0, 0.1*lc};
    topPoints[] += {4+i+1};
    
    // bottom without first and last points
    If (i > 0 && i < N)
        Point(4 + N + 1 + i) = {x0+x, -y, 0, 0.1*lc};
        botPoints[] += {4+N+1+i};
    EndIf
EndFor

botPoints[] += {4+N+1};

// rectangle lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Spline(5) = topPoints[];
Spline(6) = botPoints[];




Line Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, -6};

Plane Surface(1) = {1, 2};


// physical groups

Physical Curve("Bottom") = {1};
Physical Curve("Right") = {2};
Physical Curve("Top") = {3};
Physical Curve("Left") = {4};

Physical Curve("Airfoil") = {5,-6};

Physical Surface("Domain") = {1};

