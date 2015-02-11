Mesh.Algorithm3D = 4; // frontal
Mesh.Optimize = 1;

Mesh.CharacteristicLengthFactor = 0.15;
Merge "U_Pipe_Tri.stl";
CreateTopology;

//Extrude with grading.We define lists for the layers and heights

//r=(L-h0)/(L-hf);
//N=Ceil((Log(hf/h0))/(Log(r)));
r=1.1;
n=10;
L=0.10;

a = L*(1-r)/(1-r^n);

cells[0] = 1;
heights[0] = a;

For i In {0:n-1}
  Printf("%g",i);
  cells[i] = 1;
  heights[i] =a * (1-r^(i+1))/(1-r);

  Printf("%g",heights[i]);
EndFor
out1[] = Extrude{Surface{0}; Layers{cells[],heights[]}; Using Index[0];Recombine;};

Printf("inward extrusion: top_surf=%g volume=%g lateral_surf=%g,%g,%g", 
       out1[0], out1[1], out1[2], out1[3], out1[4]);
b1[] = Boundary{ Surface{out1[2]}; };
b2[] = Boundary{ Surface{out1[3]}; };
b3[] = Boundary{ Surface{out1[4]}; };
Printf("interior curves of lateral surfaces = %g %g %g",
       b1[2], b2[2], b3[2]);

//create inlet faces
Line Loop(200)={b1[2]};
Plane Surface(201)={200};

Line Loop(300)={b2[2]};
Plane Surface(301)={300};

//create outlet face
Line Loop(400)={b3[2]};
Plane Surface(401)={400};

//create inside volume
Surface Loop(999)={20,201,301,401};
Volume(1000) = {999};



//save only physicals
Physical Surface("inlet")={out1[2],201, out1[3], 301};
Physical Surface("outlet")={out1[4], 401};
Physical Surface("wall")={1};
Physical Volume("inside")={1000,1};
