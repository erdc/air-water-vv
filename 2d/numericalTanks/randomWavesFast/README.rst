Fast Random wave  generation
====================================

Description
----------
The test case setup is described in https://github.com/erdc-cm/air-water-vv/tree/adimako/merge/2d/numericalTanks/randomWaves.
In this particular case, the generation methodology is optimised to allow fast generation of long non repeating random waves sequences using processing with spectral windows. More details on the methodology can be found in Dimakopoulos et al. (2017).

Tests (to be added)
-----

The python test file named ``test_randomWaves.py`` is made up of two tests:

* The first test checks that the run is completed successfully.
* The second test is to validate the generation of random waves sequence at the exit of the generation zone by comparing it to the theoretical solution. 

One can run this test file typing ``py.test --boxed test_randomWaves.py``.

References
----------

- Dimakopoulos A.S., de Lataillade T and Kees C.E. Random waves in CFD models: a methodology for fast generation of long sequences, using first and second order theory, in preparation
  







