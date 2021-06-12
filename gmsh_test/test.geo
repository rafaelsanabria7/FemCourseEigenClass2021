// Gmsh project created on Sun May 30 19:19:59 2021
//+
Point(1) = {0, 0, 0, 1.0/10};
Point(2) = {1., 0, 0, 1.0/10};
Point(3) = {1, 1, 0, 1.0/10};
Point(4) = {0, 1., 0, 1.0/10};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Support_line", 5) = {1};
//+
Physical Curve("Load_line", 6) = {3};
//+
Physical Surface("Elements", 7) = {1};
