Wave sloshing
==================

Description
-----------

.. class:: center

In this case, we use Proteus to model nonlinear free-surface motion due to sloshing in a rectangular container. 
The computational domain is a 2D rectangular box with dimensions 0.1m x 0.1m and the mean water level is in the middle of the
box. The output of Proteus can be compared with the analytical solution found in Tadjbakhsh and Keller, 1960. The initial conditions of the simulation are shown in the following figure, where "a" is the amplitude of the sloshing wave.

.. figure:: ./wavesloshing.png
   :width: 100%
   :align: center

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

