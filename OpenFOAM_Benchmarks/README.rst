=====================================================
OpenFOAM Benchmarks
=====================================================

https://github.com/erdc-cm/air-water-vv/OpenFOAM_Benchmarks/

Benchmark descritpion
----------------------------

dambreakColagrossi, see  https://github.com/erdc-cm/air-water-vv/2d/dambreak_Colagrossi/README.rst
dambreakWithObstacle, see http://foam.sourceforge.net/docs/Guides-a4/UserGuide.pdf (section 2.3) and https://github.com/erdc-cm/air-water-vv/2d/dambreak_Ubbink/README.rst
lidDrivenCavity, see http://foam.sourceforge.net/docs/Guides-a4/UserGuide.pdf (Section 2.1)
dambreakGomez, see https://github.com/erdc-cm/air-water-vv/blob/master/3d/dambreak_Gomez/README.rst

Running OpenFOAM cases
------------------------
In general, how to run OpenFOAM cases is detailed in OpenFOAM's user guide (http://foam.sourceforge.net/docs/Guides-a4/UserGuide.pdf), Section 2. Regarding this test set, see recommenations below:

-For saving space, the mesh needs to be generated in all cases by running blockMesh
-The dambreakColagrossi case has a series of probes placed at the left (downstream) wall
-The dambreakGomez case has a grid of probes placed at the front and the back face of the obstacle
-The dambreakGomez case mesh is about ~6 million cells. In order to set a coarser mesh resolution, the constant/polymesh/blockMeshDict file needs to be edited before running blockMesh.

Contributors
------------
- Chris Kees, Aggelos Dimakopoulos, Eleni Zve, Matthew Farthing, Aron Ahmadia, Ido Akkerman, Matt Malej, Roham Bakhtyar


