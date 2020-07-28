SetFactory("OpenCASCADE");


x1 = 2*Pi*256/1000;
x2 = 2*Pi*256/1000+2*Pi;
n =1000;
lc1 = 0.085;
lc2 = 0.085;
lc3 = 0.085;
lc4 = 0.085;

//First interface Points and Lines

For i In {0:n-1}

  x = x1 + (x2 - x1) * i/n;

  r = 0.7512+0.0124*Cos(x)+0.0414*Cos(2*x)+0.0120*Cos(3*x)+0.0111*Cos(4*x)-0.0061*Cos(5*x)-0.0008*Cos(6*x)-0.00001*Cos(7*x)+0.0033*Cos(8*x)+ 0.0217*Sin(x)-0.0197*Sin(2*x)-0.0027*Sin(3*x)-0.0084*Sin(4*x)-0.0135*Sin(5*x)-0.0034*Sin(6*x)-0.0001*Sin(7*x)-0.0048*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), 0, lc1}; //if you want the inner boundary mesh to be more refined change lc1 to a smaller number

EndFor

s1 = newreg;
Spline(s1) = {pList1[{0:n-1}]};

//Second interface Points and Lines
x1 = 2*Pi*255/1000;
x2 = 2*Pi*255/1000+2*Pi;
For i In {n:2*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = 1.0284+0.0362*Cos(x)-0.0557*Cos(2*x)+0.0190*Cos(3*x)-0.0185*Cos(4*x)+0.0038*Cos(5*x)+0.0045*Cos(6*x)-0.0002*Cos(7*x)+0.0006*Cos(8*x)-0.1984*Sin(x)-0.0173*Sin(2*x)-0.0008*Sin(3*x)-0.0133*Sin (4*x)-0.0026*Sin(5*x)-0.0066*Sin(6*x)-0.0029*Sin(7*x)+0.0064*Sin(8*x);


  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), 0, lc2};

EndFor

s2 = newreg;
Spline(s2) = {pList1[{n:2*n-1}]};


//Third Interface Points and Lines
x1 = 2*Pi*254/1000;
x2 = 2*Pi*254/1000+2*Pi;
For i In {2*n:3*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = 1.2584+0.0362*Cos(x)-0.0557*Cos(2*x)+0.0190*Cos(3*x)-0.0185*Cos(4*x)+0.0038*Cos(5*x)+0.0045*Cos(6*x)-0.0002*Cos(7*x)+0.0006*Cos(8*x)-0.1984*Sin(x)-0.0173*Sin(2*x)-0.0008*Sin(3*x)-0.0133*Sin(4*x)-0.0026*Sin(5*x)-0.0066*Sin(6*x)-0.0029*Sin(7*x)+0.0064*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), 0, lc3};

EndFor


s3 = newreg;
Spline(s3) = {pList1[{2*n:3*n-1}]};



//Fourth interface Points and Lines
x1 = 2*Pi*254/1000;
x2 = 2*Pi*254/1000+2*Pi;
For i In {3*n:4*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = 1.3584+0.0362*Cos(x)-0.0557*Cos(2*x)+0.0190*Cos(3*x)-0.0185*Cos(4*x)+0.0038*Cos(5*x)+0.0045*Cos(6*x)-0.0002*Cos(7*x)+0.0006*Cos(8*x)-0.1984*Sin(x)-0.0173*Sin(2*x)-0.0008*Sin(3*x)-0.0133*Sin(4*x)-0.0026*Sin(5*x)-0.0066*Sin(6*x)-0.0029*Sin(7*x)+0.0064*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), 0, lc4};

EndFor


s4 = newreg;
Spline(s4) = {pList1[{3*n:4*n-1}]};

//The following lines should be written by hand. For extracting the labels you have to check them on GMSH. For example the number {2} for the first curveloop was 
//exracted by hovering the mouse pointer over the inner interface.



//+
Line(5) = {1, 1001};
//+
Line(6) = {1001, 2001};
//+
Line(7) = {2001, 3001};
//+
Line(8) = {1000, 2000};
//+
Line(9) = {2000, 3000};
//+
Line(10) = {3000, 4000};
//+
Curve Loop(1) = {1, 8, -2, -5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 9, -3, -6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, 10, -4, -7};
//+
Plane Surface(3) = {3};
//+
Extrude {0, 0, 1} {
  Curve{1}; Surface{1}; Curve{5}; Curve{8}; Curve{2}; Surface{2}; Curve{6}; Curve{9}; Curve{3}; Surface{3}; Curve{7}; Curve{10}; Curve{4}; 
}


//+
Physical Surface(4) = {8,12,16};
//+
//Physical Surface(5) = {12};
//+
//Physical Surface(6) = {16};
//+
Physical Surface(5) = {1,2,3};
//+
//Physical Surface(8) = {2};
//+
//Physical Surface(9) = {3};
//+
Physical Surface(6) = {7,11,15};
//+
//Physical Surface(11) = {11};
//+
//Physical Surface(12) = {15};
//+
Physical Surface(7) = {5,9,13};
//+
//Physical Surface(14) = {9};
//+
//Physical Surface(15) = {13};
//+
Physical Surface(8) = {4};
//+
Physical Surface(9) = {6};
//+
Physical Surface(10) = {10};
//+
Physical Surface(11) = {14};
//+
Physical Volume(12) = {1};
//+
Physical Volume(13) = {2};
//+
Physical Volume(14) = {3};
//+

