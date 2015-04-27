Wave diffraction - Penny and Price (1952)
==================================================================

Wave diffraction presents a boundary behaviour of waves associated with the bending of their path. Wave diffraction involves a change in direction of the wave propagation as they pass through an opening or around a barrier. The waves are, then, seen to pass around the barrier into regions behind it and disturb the water behind it. 
The present benchmark case represents the case of wave diffraction caused by a semi-infinite breakwater. Analytical solutions for this type of diffraction are presented in Penny and Price (1952) and can be used for comparisons of the corresponding results of PROTEUS. 
The case is simulated using a 3D rectangular numerical domain with height of 1.5 m, a length, L, of 25.0 m and width, b, of 30.0 m. The breakwater is represented as a vertical thin solid wall with length, b1, equal to 10.0 m. The mean water depth within the domain is equal to 1.0 m. within the left boundary and the generation/absorption zone, non-linear waves are generated with a height of 0.10 m, a period of 1.94 s and a length of 5.0 m. Also, the waves are absorbed within the implemented absorption zone. For the determination of the flux parameters, the Fenton Fourier Transformation theory is used (Fenton, 1988). Atmospheric conditions have been assigned to the top boundary of the domain and the bottom boundary, the surfaces of the breakwater, the front and the back boundaries act as free-slip walls.

A sketch of the plan of the 3D domain of the present case is given in the following figure.

.. figure:: ./Diffraction_PP_plan.bmp

where, 

* Lga=1 wavelength (=5.0m), the length of the Generation/absorption zone 
* L1=L2=1 wavelength (=5.0m)
* La=2 wavelengths (=10.0m), the length of the Absorption zone in the x direction
* ba=2 wavelengths (=10.0m), the length of the Absorption zone in the y direction
* b1=b2=2 wavelengths (=10.0m)
       
This test case demonstrates the ability of PROTEUS to simulate the wave diffraction caused by a semi-infinite breakwater as well as the wave absorption within a 3D domain configuration.

References
--------------------------------

- Penny WG and Price AT (1952) The diffraction theory of sea waves and the shelter afforded by breakwaters. Part I in \u201cSome gravity wave problems in the motion of perfect liquids\u201d by Martin JC, Moyce WJ, Penny WG, Price AT and Thornhill CK, Philosophy Transaction of the Royal Society in London, A244, 236\u2013253.

- Fenton JD (1988) The numerical solution of steady water wave problems. Computer and Geosciences, 14(3), 357-368.


