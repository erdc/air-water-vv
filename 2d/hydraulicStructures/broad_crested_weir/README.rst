Broad crested weir
==================

Description
-----------

A broad-crested weir is a standard hydraulic structure used as
discharge measuring device and flow control device in open
channel. This type of weir can be described as a simple solid
rectangular channel obstruction spanned over the whole width of the
channel. Due to the sharp edge of the upper left side corner of the
weir, flow separation occurs at this location.

As the flow propagates mainly in 2 directions, a 2D computational
domain was used for the simulation.  The height of the domain is equal
to :math:`1.0 m` and the length to :math:`3.5 m`.
The weir has a height of :math:`P=0.401 m` and a width of
:math:`b=0.5 m`.  A uniform velocity distribution from bottom to
free water level is imposed in the left wall boundary condition. The
top of the domain is left open and the right wall of the domain allows
the flow to leave the domain. A free-slip condition is set at the 
bottom of the domain. There is an air ventilation on the right side of 
the weir to maintain atmospheric pressure below the nappe. In the following figure, a simple sketch of the structure is 
presented showing the main parameters.

.. figure:: ./BroadWeir.bmp
   :width: 100%
   :align: center

where, :math:`u_0` is the approach velocity, :math:`H` is the upstream
potential head, :math:`hv` is the upstream velocity head, :math:`Ht =
H + hv` is the upstream total head, :math:`d` is the flow depth over
the weir.

This case tests the ability of Proteus to simulate the free-surface
evolution and the flow separation. The results of the simulations can
be compared with the data in the following references.  For more
details, see runfiles or references.


Running the test case
-----

To run the test case type:

```
parun broad_crested_weir.py --TwoPhaseFlow -v -D result_folder
```

Geometry and setup options can be modified by the commandline, using for example:

```
parun broad_crested_weir.py --TwoPhaseFlow -v -D result_folder -C "obstacle_dim=(0.5,0.4)"
```

To run in parallel (example with mpirun and 12 processors):

```
mpirun -np 12 parun broad_crested_weir.py --TwoPhaseFlow -v -D result_folder -C "obstacle_dim=(0.5,0.4)"
```


To see guidance on parun options, you can type  

```
parun -h
```


References
----------

- Fritz HM and Hager WH (1998) Hydraulics of embankment weirs. Journal
  of Hydraulic Engineer 124(9), 963–971.

- Hager WH and Schwalt M (1994). Broad-crested weir. Journal of
  Irrigation and Drainage 120(1), 13–26.

