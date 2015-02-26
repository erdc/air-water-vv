//The isotropic mesh size can be prescribed by defining
//a mesh characteristic length factor
//Mesh.CharacteristicLengthFactor=1;
Merge "U_pipe.step";
lc=0.01;

Characteristic Length {1}=lc;
For num In {1:20}
		Characteristic Length {num+1}=lc;
EndFor

Mesh.CharacteristicLengthMax=lc;
//In case the triangulation is not a closed surface the topology
//of the mesh has to be created
CreateTopology;

//Characteristic Length {1,2,3,4}=lc;
//Characteristic Length {5,6,7}=lc;
Mesh.Smoothing=5;
//Mesh.Algorithm=5;	// 2D Mesh Algorithm
Mesh.Algorithm3D=1; // 3D Mesh Algorithm
Mesh.SaveElementTagType = 2;
Mesh.QualityType = 1;
Mesh.Optimize=1;

//Field[7] = BoundaryLayer;
//Field[7].FacesList = {1,3};
//Field[7].Quads = 1;
//Field[7].hfar = 1;		// element size far from wall
//Field[7].hwall_n = ;	// mesh size normal to the wall


Physical Surface(3) = {1, 3, 5};
Physical Surface(1) = {2};
Physical Surface(2) = {4};
Physical Volume(1) = {1};

Geometry.Tolerance = 1e-2;
Geometry.OCCSewFaces = 1;
