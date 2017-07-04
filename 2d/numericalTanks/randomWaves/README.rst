Random waves generation/absorption
====================================

Random waves typically consist of non-repeatable wave sequences, so each individual wave has different characteristics. Random waves in nature usually obey spectral distributions e.g. JONSWAP and Pierson-Moskowitz spectra and their statistical properties can be predicted. An example of a spectral distribution is shown in the figure below.

.. figure:: ./Spectrum.PNG
   :width: 100%
   :align: center


The python test file named ``test_nonlinearWaves.py`` is made up of three tests:

* The first test checks that the run is completed successfully.
* The second test is to validate the results comparing them to the theory. For this case we will compare the numerical and theoretical wave height in the middle of the tank.
* The third test evaluates wave reflection and compares to a threshold. The calculation of reflection is performed by applying Isaacson's 3rd method (Isaacson 1991) to the primary harmonic of the signal.

One can run this test file typing ``py.test --boxed test_nonlinearWaves.py``.

References
----------

- Fenton JD (1988) The numerical solution of steady water wave 
  problems, Comp and Geosc, 14(3), 357-368
  
- Lé Méhauté, B., (1976). “Introduction to Hydrodynamics and water waves”, Springer-Verlag, New York.

- Isaacson (1991), Measurement of regular wave reflection, Journal of Waterway Port Coastal and Ocean Engineering 117(6), 553-569





