backres = DefineNumber[ 1000, Name "Parameters/backres" ];
frontres = DefineNumber[ 1000, Name "Parameters/frontres" ];
Point(1) = {0, 0, 0, backres};
Point(2) = {85000, 0, 0, backres};
Point(3) = {90000, 0, 0, frontres};
Point(4) = {105000, 0, 0, frontres};
Point(5) = {110000, 0, 0, backres};
Point(6) = {480000, 0, 0, backres};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Extrude {0, 20, 0} {
  Line{1}; Line{2}; Line{3}; Line{4}; Line{5}; Layers{10}; Recombine;
}


Physical Line(1) = {1, 2, 3, 4, 5};
Physical Line(2) = {6, 10, 14, 18, 22};
Physical Line(3) = {7};
Physical Line(4) = {24};
Physical Surface(1) = {9, 13, 17, 21, 25};
