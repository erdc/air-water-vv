Dambreak flow - Collagrosi and Landrini (2003)
==============================================

Description
-----------
The problem comprises a 2D tank with dimensions  3.22m x 1.8m  (width x height). 
The boundaries are considered closed, specifically as free slip walls, while the top boundary [y+] is set as an atmposphere.
The purpose of this case initially published by Collagrosi and Landrini (2003) is to examine a two phase flow in a 2D tank, where the 2 fluids are separated by sharp interfaces.
The initial condition of the fluid is a water column with dimensions  1.20m x 0.60m  (width x height). By the time the simulation starts, the column collapses under the action of gravity and impacts to a wall.

In the following figure, the geometry of the dambreak case is illustrated.

.. figure:: ./dambreakColagrossi.bmp
   :width: 100%
   :align: center

This case tests the ability of Proteus to simulate the free-surface
evolution and forces / pressures on structures, according to data that
are available in the following references.  For more details, see
runfiles or references.

Running the test case
-----

To run the test case type:

```
parun dambreak.py--TwoPhaseFlow  -v -D result_folder
```

Dambreak and tank properties can be modified by the commandline, using for example:

```
parun dambreak.py --TwoPhaseFlow  -v -D result_folder -C "mwl=0.5"
```

To run in parallel (example with mpirun and 12 processors):

```
mpirun -np 12 parun -f dambreak.py --TwoPhaseFlow  -v -D result_folder -C "mwl=0.5"
```


To see guidance on parun options, you can type  

```
parun -h
```


References
----------

- Colagrossi A and Landrini M (2003) Numerical simulation of
  interfacial flows by smoothed particle hydrodynamics, Journal of
  Computational Physics,191,448-475.

- Martin, J. C. & Moyce, W. J., (1952) Part IV. An Experimental Study
  of the Collapse of Liquid Columns on a Rigid Horizontal Plane
  Phil. Trans. R. Soc. Lond. A 244 (882) 312-324.

- Zhou, Z. Q., De Kat, J. O. and Buchner, B. (1999) A nonlinear 3-D
  approach to simulate green water dynamics on deck in: J. Piquet
  (Ed.), Proc. 7th Int. Conf. Num. Ship Hydrod., Nantes, 5.11, 15.
