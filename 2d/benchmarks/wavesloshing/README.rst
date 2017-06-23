Wave sloshing
==================

Proteus is initialised with a free-surface profile that will excite
the first asymmetric sloshing mode. The free-surface slope is small,
therefore, we can consider that the flow evolution will obey to the
linear sloshing theory. The output of Proteus can be compared with the
analytical solution found in Standing surface waves of finite amplitude
(Tadjbakhsh and Keller, 1960) including only a symmetric
mode. The computational domain is a 2D rectangular box with dimensions
0.1m x 0.1m and the mean level of the water is in the middle of the
box. The initail conditions of the simulation are shown in the
following figure,

.. figure:: ./wavesloshing_scheme.png
   :width: 100%
   :align: center

where, a is the amplitude of the sloshing wave.

This case tests the ability of PROTEUS to simulate the free-surface
evolution. For more details, see runfiles or references.

The python test file named ``test_wavesloshing.py`` is made up of 
two tests:

* The first one is to know if the case can run.
* The second test is to validate the results comparing them to the theory. For this case we will compare the numerical and theoretical position of the interface at the left boundary.
One can run this test file typing ``py.test --boxed test_wavesloshing.py``.

References
----------

- Tadjbakhsh, I., & Keller, J. (1960). Standing surface waves of finite amplitude. Journal of Fluid Mechanics, 8(3), 442-451. doi:10.1017/S0022112060000724

