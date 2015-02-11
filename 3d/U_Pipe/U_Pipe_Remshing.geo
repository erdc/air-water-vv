//When remeshed, a compound surface will be reparametrized with either harmonic or conformal maps.You can specify this as follows at the beginning of the geo file: 
Mesh.RemeshParametrization=0; //0=harmonic_circle, 1=conformal_spectral, 2=rbf, 3=harmonic_plane, 4=convex_circle, 5=convex_plane, 6=harmonic square,7=conformal_fe 
Mesh.Optimize = 1;
Mesh.RemeshAlgorithm=1; //(0) nosplit (1) automatic (2) split only with metis

// Those algorithms are implemented in Gmsh. In order to remesh stl (also possible for .msh, .step,..) files, you first have
// to merge your stl file using the Merge command and then define compound surfaces
Mesh.CharacteristicLengthFactor = 0.1; //The isotropic mesh size can be prescribed by defining a mesh characteristic length factor.
Merge "U_Pipe.stl";
CreateTopology; //In case the triangulation is not a closed surface the topology of the mesh has to be created.

// A compound surface is defined as a set of several elementary discrete surfaces (in the case of stl triangulations, 
// there is only one discrete triangulation that has tag '1'). 
Compound Surface(1) = {1:1000}; // auto-detect boundary
Physical Surface(101)={1};