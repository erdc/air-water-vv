Crump weir
==========

A crump weir is a typical hydraulic structure mainly used as a
measurement device in open-channels as well as for flow control. It is
a triangular channel obstruction consisting of a 1:2 upstream and a
1:5 downstream slope.The slope upstream is designed to avoid material
sedimentation while the downstream side is designed to stabilize the
flow. The weir spans the width of the channel and obstructs the
natural flow of the water.The following two flow conditions are
tested:

* Modular flow conditions, where the flow regime is not affected by
  the downstream conditions. In this condition there is a direct
  relation between the flow rate and the upstream head. The modular
  flow condition is verified when the ratio between the downstream and
  the upstream total head is less or equal to 0.75.

* Non-modular flow conditions, where the flow regime is affected by
  the downstream conditions. This happens when the weir operates in
  submerged conditions.  The flow rate in this case is related to the
  total head upstream and downstream of the weir and therefore both
  head measurements are required.

The computational domain is a 2D rectangular box with height equal to
2.1m the length might vary from case to case.The weir has a height of
:math:`P=0.5\mbox{m}` and a width of :math:`b=3.5\mbox{m}`.  A uniform
velocity distribution from bottom to free water level is imposed in
the left wall boundary condition.  The top of the domain is left
open. The right wall of the domain allows the flow to leave the domain
for the case of the modular flow, when for the case of the non-modular
flow an outflow velocity and downstream water level are imposed on the
right boundary.

In the following figure a simple sketch of the structure is presented.

.. figure:: ./CrumpWeirBSISO.jpg
   :width: 100%
   :align: center

This case tests the ability of PROTEUS to simulate the free-surface
evolution and the results of the simulations can be compared with the
data in the following references.  For more details, see runfiles or
references.

The python test file named ``test_crump_weir.py`` is made up of 
two tests:
* The first one is to know if the case can run.
* The second test is to validate the results comparing them to the theory. 
For this case we will compare the numerical and theoretical discharge 
over the crest of the weir.
One can run this test file typing ``py.test --boxed test_crump_weir.py``.

References
----------

- British Standards Institution (2008) BS ISO 4360:2008: “Hydrometry -
  Open channel flow measurement using triangular profile
  weirs”. London, BSI. (Withdrawn (ISO 4360:1984 is a current
  alternative).
