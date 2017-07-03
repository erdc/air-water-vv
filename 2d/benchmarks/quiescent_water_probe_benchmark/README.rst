Quiescent water
==================

Description
-------------
The computational domain consists of a 2D rectangular box, initialised 
with a quiescent free surface, so that the computational domain includes
both the air and the water phase.  The flow variables (pressure and 
velocities) are not expected to vary in space and time, as the 
simulation progresses in time. The computational domain is a 2D 
rectangular box with dimensions 3.22m x 1.8m and the level of water is 
at 0.6m. The initial conditions of the simulation are shown in the
following figure.

.. figure:: ./quiescent_water_scheme.png
   :width: 100%
   :align: center
   
Tests
-------
The python test file named ``test_quiescent_water.py`` is made up of 
two tests:

* The first test is to check that the run is successfully completed.
* The second test is to validate the readings of the probes by comparing them to theoretical values. For this case we will compare the numerical and theoretical pressure and water level.
One can run this test file typing ``py.test --boxed test_quiescent_water.py``.
