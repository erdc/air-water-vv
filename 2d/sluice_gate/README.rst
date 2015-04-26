Sluice gate
===========

An under-shot sluice gate is a hydraulic structure mainly used as a
flow control device in open-channels, where the discharge downstream
can be related to the upstream head and therefore controlled.  It is
basically a vertical obstruction spanning the whole width of the
channel with an opening across the bottom.  The computational domain
is a 2D rectangular box with height equal to 3.75m and the length
equal to 4.5m The sluice gate has an opening of :math:`P=0.25\mbox{m}`
and a small width of :math:`b=0.01\mbox{m}`. A uniform velocity
distribution from bottom to free water level is imposed in the left
wall boundary condition. The top of the domain is left open and the
right wall of the domain allows the flow to leave the domain.  In the
following figure a simple sketch of the structure is presented showing
the main parameters.

.. figure:: ./SluiceGate.bmp
   :width: 100%
   :align: center

where, :math:`u_0` is the approach velocity, :math:`u_1` is the
velocity downstream of the gate, d1 is the upstream water level, P is
the opening of the gate, d1 is the water depth downstream of the gate.

This case tests the ability of PROTEUS to simulate the free-surface
evolution.  For more details, see runfiles or references.


References
----------

- White F.M. 1999. “Fluid Mechanics”. Fourth Edition, McGraw-Hill
  series in mechanical engineering.

