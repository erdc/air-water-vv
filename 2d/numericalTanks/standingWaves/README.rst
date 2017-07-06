Standing Waves
====================================

Description
----------

Standing waves are formed when regular waves interact with a fully reflective wall. In this case, a distinctive patterns emerges and the wave amplitude becomes a sinusioidal function in space. The reflected waves are fully synconised with the incident ones at an even number of half wavelenghts from the wall, whilst cancelling each other at an odd number of wavelength, this creating the distinctive wave patterns of nodes and antinodes, where the amplitude is double and zero, respectively.

In this case, the numerical flume described in https://github.com/erdc-cm/air-water-vv/tree/adimako/merge/2d/numericalTanks/linearWaves is used to demonstrate the capability of Proteus of modelling standing wave patterns and in particular of absorbing the reflected waves at the generation / absorption zone.

Tests
-----

The python test file named ``test_standingWaves.py`` is made up of three tests:

* The first test checks that the run is completed successfully.
* The second test is to validate the model against the analytical solution of the free-surface elevation at the reflective wall

One can run this test file typing ``py.test --boxed test_standingWaves.py``.







