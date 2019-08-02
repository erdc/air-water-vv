Wave sloshing
==================

Description
-----------

In this case, we use Proteus to model nonlinear free-surface motion due to
sloshing in a rectangular container. The output of Proteus is compared with the
analytical solution found in Tadjbakhsh and Keller, 1960. The computational domain is a 2D rectangular box with dimensions
0.1m x 0.1m and the mean level of the water is in the middle of the
box. The initail conditions of the simulation are shown in the
following figure.

.. figure:: ./wavesloshing.png
   :width: 100%
   :align: center

where, a is the amplitude of the sloshing wave.

This case tests the ability of PROTEUS to simulate the free-surface
evolution. For more details, see runfiles or references.

Test case
-----

The test case comprises a simple rectangular tank with generation zone at the left side ('x-') and absoprtion zone at the right side ('x+'). To run the test case type:

```
parun --TwoPhaseFlow -f wavesloshing.py -v -D result_folder
```

Wave properties can be modified by the commandline, using for example:

```
parun --TwoPhaseFlow -f wavesloshing.py -v -D result_folder -C "T=10."
```

To run in parallel (example with mpirun and 12 processors):

```
mpirun -np 12 parun --TwoPhaseFlow -f wavesloshing.py -v -D result_folder -C "T=10."
```


To see guidance on parun options, you can type  

```
parun -h
```

References
----------

- Tadjbakhsh, I., & Keller, J. (1960). Standing surface waves of finite amplitude. Journal of Fluid Mechanics, 8(3), 442-451. doi:10.1017/S0022112060000724

