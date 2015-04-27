Wave breaking - Ting and Kirby benchmark
=============

Sea waves transform when propagating from the offshore to the near
shore environment and as they interact with coastal structures. In a
2D flume, the main factors contributing to wave transformation are
shoaling, reflection, transmission and breaking. The present case aims
to reproduce numerically the laboratory experiments on wave breaking
of Ting and Kirby (1994).  Ting and Kirby (1994) conducted experiments
on a 1 to 35 constant-slope bathymetry on 2 different types of
breaking i.e. spilling and plunging breaking of cnoidal waves.

The present problem consists of a 2D numerical flume with a sloping
bottom of 1:35 slope. Upstream the slope the domain has a height, H,
of 1.0 m, where the mean water depth, h, is equal to 0.4 m. The total
length of the domain, L, is 40m. Two cases are simulated with regard
to the type of breaking i.e. (a) spilling and (b) plunging
breaking. In both cases, at the left boundary cnoidal waves are
generated. The wave heights are 0.125 m and 0.128 m, for cases (a)
and (b) respectively, when the wave period in the case (a) is 2.0s and
in case (b) is 5.0s. For the generation of the waves the Fenton's
method (Fenton, 1988) have been used for identifying the input flux
parameters within the left boundary. The bottom and right boundaries
act as a free-slip walls, when in the top boundary have been assigned
atmospheric conditions. A sketch of the domain is given in the
following figure.

.. figure:: ./Breaking.bmp
   :width: 100%
   :align: center

where L1=15.0m and ht=0.38m

This test case demonstrates the ability of PROTEUS to simulate the
shoaling process of two different types of wave breaking over a
constant slope bathymetry.

References
----------

- Ting FCK and Kirby JT (1994) Observation of Undertow and Turbulence
  in a Laboratory Surf Zone. Coastal Engineering, 24, 177-204.

- Fenton JD (1988) The numerical solution of steady water wave
  problems. Computer and Geosciences, 14(3), 357-368.


