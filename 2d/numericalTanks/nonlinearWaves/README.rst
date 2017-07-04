Nonlinear wave generation/absorption
====================================

Plane regular nonlinear waves are mild and high steepness waves that 
propagate in a single direction, in uniform wave fronts.  The wave 
profile deviates from the sinusoidal shape, and it typically exhibits 
high and sharp wave crests and low and flat wave troughs.  It is not 
always accurate to calculate the wave celerity from the linear 
dispersion theory, especially for highly nonlinear waves. 
Fenton (1988) proposes a method for calculating the nonlinear wave 
properties and profile, which is adopted for the generation of 
nonlinear waves within Proteus. 

In terms of classification, linear waves are in the right top area of the following diagram (Lé Méhauté 1976).    where, the vertical axis corresponds to the no dimensional wave height 
and the horizontal to the no dimensional water depth. The term gT\ 
:sup:`2`\ is proportional to the wavelength in deep water and the dot 
named A corresponds to the tested case which is described below.


.. figure:: ./Mehaute_nonlinear_waves_01.png
   :width: 100%
   :align: center

  The 
The numerical wave flume represents the geometry used in Higuera et al 2013 for their numerical tests. The file nonlinearTest
This case tests demonstrates the ability of PROTEUS to simulate the 
generation and propagation of non-linear waves as well as their 
absorption. 

The python test file named ``test_nonlinearWaves.py`` is made up of three tests:

* The first one is to know if the case can run.
* The second test is to validate the results comparing them to the theory. For this case we will compare the numerical and theoretical wave height in the middle of the tank.
* The third one is to test the reflection. 
One can run this test file typing ``py.test --boxed test_nonlinearWaves.py``.

References
----------

- Fenton JD (1988) The numerical solution of steady water wave 
  problems, Comp and Geosc, 14(3), 357-368.







