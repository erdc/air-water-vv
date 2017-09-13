Flat plate turbulent flow – Wall function benchmark
==============================================

Description
-----------
The problem comprises a 0.40m x 4.00m (height x length)tank with flat plate either at the bottom and the top.
The water flows through the duct and interacts with the solid surfaces generating turbulence in the near-wall region. 

.. figure:: ./flatPlateBenchmark.bmp
   :width: 100%
   :align: center

At the very proximity of the wall, in the viscous sub-layer, viscous contribution to the shear stress is significant and in general Reynolds stresses are negligible when compared with it.
Turbulent effects become gradually more important in the inner region moving away from the solid surface. Velocity field assumes a logarithmic profile and the viscous effects can be considered negligible at some point (log-law sublayer).

Simplified approaches to this problem are adopted and the near-wall sublayer is not resolved. Empirical formulations of the problem based on experimental data observations are used. The law of the wall for the inner part of a turbulent shear flow over a solid surface includes a simple analytic function for the mean velocity distribution, the logarithmic law.

This case tests the ability of PROTEUS to simulate the generation of turbulence at solid surfaces coupled with velocity assuming logarithmic profile in the inner region.
For more details, see references.


References
----------

- Pope S.B., Turbulent Flows. Wall flows, 264–298. Reynolds-stress and related models, 442-444.

- Schlichting H., Boundary Layer Theory. Turbulent flow through pipes, 596-623.

