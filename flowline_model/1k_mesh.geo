backres = DefineNumber[ 1000, Name "Parameters/backres" ]; 
frontres = DefineNumber[ 1000, Name "Parameters/frontres" ]; 
Point(1) = {0.0, -3000.0, 0, backres}; 
Point(2) = {480000.0, -3000.0, 0, backres}; 
Line(1) = {1, 2};
Extrude {0, 3000, 0} {
  Line{1};Layers{5};Recombine;
}
Extrude {0, 20, 0} {
  Line{2};Layers{10};Recombine;
}

Physical Surface(1) = {9};
Physical Surface(3) = {5};
Physical Line(1) = {2};
Physical Line(2) = {6};
Physical Line(3) = {7};
Physical Line(4) = {8};
Physical Line(5) = {1};
Physical Line(6) = {3};
Physical Line(7) = {4};

