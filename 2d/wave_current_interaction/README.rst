Wave/current interaction
========================

The presence of a current in the wave field causes the waves to
transform. The primary effect of a current, with velocity u, to a wave
is the change the wavelength due to a Doppler shift. If the wave and
current travel in the same direction (following current), then the
wavelength is increased over the one that waves would possess in still
water. If waves and current travel in opposite direction (opposing
current) then the wavelength is decreased. Wave current interaction
with the waves also affects the wave height due to wave energy
conservation and the wave height decreases for following currents and
increases for opposing currents. Given the wave period, modified
height, depth and current velocity, the Fenton Fourier Transform
theory can be used to calculate wavelength, particle velocities and
pressure and these variables can be imposed as a boundary conditions
in a numerical wave flume.

A sketch of the interaction of waves with a current is given in the
following figure.


.. figure:: ./  WaveCurrent.bmp
   :width: 100%
   :align: center

where, lambda0, lambda1, lambda2 are the wave lengths in the case
without a current, an opposing current and a following current,
respectively and the H0, H1 ,H2 are the corresponding wave heights.
Also, the relationship between the wave lengths and the wave heights
is the following: lambda1 < lambda0 < lambda2 and H2 < H0 < H1.

The present benchmark case aims to simulate 2 different cases of
wave/current interaction, where the modification of a linear regular
wave by a uniform steady current is modelled. The first case consists
of a regular wave with a following current of 1 m/s and the second
case to an opposing current of 0.5 m/s to the same regular wave as in
the first case. The numerical domain of both cases consists of a 2D
rectangular numerical flume with height of 1.50 m and a length depending 
on the wavelength. The mean water depth is equal to 1.0 m. The flux parameters 
of the modified waves due to a uniform current are imposed at the left boundary 
of the domain and the flux parameters are defined using Fenton's method (Fenton, 1988). 
Atmospheric conditions have been assigned to the top boundary of the domain, 
the bottom boundary acts as a free-slip wall.

This test case demonstrates the ability of PROTEUS to simulate the
wave/current interaction in a 2D configuration.

References
----------

- Brevik I and Aas B (1980) Flume experiments on waves and
  currents II. smooth bed. Coastal Engineering,3, 149-177.

- Fenton JD (1988) The numerical solution of steady water wave
  problems, Comp and Geosc, 14(3), 357-368.

