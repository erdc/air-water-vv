Dambreak flow - Collagrosi and Landrini (2003)
==============================================

Description
-----------
The problem comprises a 0.60m x 1.20m (height x width) column of
water in a 1.8 m high and 3.22 m wide container. The column collapses under the action of gravity
and impacts to a wall. The top of the domain is
left open, when the rest of the boundary patches act as free slip walls.
In the following figure, a sketch of the dambreak initial conditions
is shown.

.. figure:: ./dambreakColagrossi.bmp
   :width: 100%
   :align: center

This case tests the ability of Proteus to simulate the free-surface
evolution and forces / pressures on structures, according to data that
are available in the following references.  For more details, see
runfiles or references.

Test case
-----

The test case comprises a simple rectangular tank with generation zone at the left side ('x-') and absoprtion zone at the right side ('x+'). To run the test case type:

```
parun --TwoPhaseFlow -f dambreak.py -v -D result_folder
```

Wave properties can be modified by the commandline, using for example:

```
parun --TwoPhaseFlow -f dambreak.py -v -C -D result_folder "Tp=2 Hs=0.2"
```

To run in parallel (example with mpirun and 12 processors):

```
mpirun -np 12 parun --TwoPhaseFlow -f dambreak.py -v -C -D result_folder "Tp=2 Hs=0.2"
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
