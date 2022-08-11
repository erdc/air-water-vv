Fast Random wave  generation
====================================

Description
----------
The test case setup is described in https://github.com/erdc/air-water-vv/tree/adimako/merge/2d/numericalTanks/randomWaves.
In this particular case, the generation methodology is optimised to allow fast generation of long non repeating random waves sequences using processing with spectral windows. More details on the methodology can be found in Dimakopoulos et al. (2017).

Tests (to be added)
-----

The python test file named ``test_randomWavesFast.py`` is made up of three tests:

* The first test is a fast test to see if the case can run for a very short simulation time. It only can be 
  launch from the air-water-vv folder typing ``py.test --boxed -v --runfast 2d/numericalTanks/randomWavesFast/test_randomWavesFast.py``.
* The second test checks that the run is completed successfully.
* The third test is to validate the generation of random waves sequence at the exit of the generation zone by comparing it to the theoretical solution. 

One can run the two last tests file typing ``py.test --boxed -v test_randomWavesFast.py`` in the randomWavesFast folder.

References
----------

- Dimakopoulos A.S., de Lataillade T and Kees C.E. Random waves in CFD models: a methodology for fast generation of long sequences, using first and second order theory, in preparation
