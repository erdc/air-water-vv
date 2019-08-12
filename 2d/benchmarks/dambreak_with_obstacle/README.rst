Dambreak flow with obstacle - Ubbink (1997)
===========================================

Description
-----------

The problem comprises a 0.292m x 1.146m (height x width) column of
water in a rectangular container (0.584 m x 0.584 m), initially at rest.  A small obstacle (0.024m x 0.048m) is
placed in the middle of the container's base. The top of the domain is left open, when the rest of
the boundary patches act as free slip walls.  In the following figure, a
sketch of the dambreak initial conditions is shown.

.. figure:: ./dam_break_wiht_obstacle.jpg
   :width: 100%
   :align: center

Under the action of gravity the water column collapses and interacts with the obstacle leading 
to a formation of a jet that impacts to a wall.

This case tests the ability of PROTEUS to simulate the free-surface
evolution and during the interaction of the dambreak front with the
obstacle.  The results of the simulations can be visually compared
with screenshots from Ubbink (1997).  For more details, see
runfiles or references.

Test case
-----

To run the test case type:

```
parun --TwoPhaseFlow -f dambreak_with_obstacle.py -v -D result_folder
```

Dambreak and tank properties can be modified by the commandline, using for example:

```
parun --TwoPhaseFlow -f dambreak_with_obstacle.py -v -D result_folder -C "tank_dim=(0.8, 0.8)"
```

To run in parallel (example with mpirun and 12 processors):

```
mpirun -np 12 parun --TwoPhaseFlow -f dambreak_with_obstacle.py -v -D result_folder -C "tank_dim=(0.8, 0.8)"
```


To see guidance on parun options, you can type  

```
parun -h
```

References
----------

- Ubbink, O. (1997), Numerical prediction of two fluid systems with
  sharp interfaces, PhD thesis, Department of Mechanical Engineering,
  Imperial College of Science, Technology & Medicine
