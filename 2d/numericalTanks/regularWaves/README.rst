Regular Wave Generation
====================================

Description
-----------

This test case is used for generating regular waves in a numerical tank, for both linear and nonlinear regimes. 

Water waves with low steepness are considered linear (:math:`H/L` < 0.1%, where :math:`H` is the wave height and :math:`L` is the wavelength). This means that the form of the waves is sinusoidal and high order terms are negligible. 
The wavelength, wave period and water depth are interdependent through the linear dispersion relation. 


Regular nonlinear waves are mild and high steepness waves that 
propagate in a single direction, in uniform wave fronts.  The wave 
profile deviates from the sinusoidal shape, and it typically exhibits 
high and sharp wave crests and low and flat wave troughs.
Fenton (1988) proposes a method for calculating the nonlinear wave 
characteristics and profile, which is adopted for the generation of 
nonlinear waves within Proteus. 

In terms of classification, linear waves are in the lower right corner of the diagram (Lé Méhauté 1976), where the vertical axis corresponds to the non dimensional wave height and the horizontal to the non dimensional water depth, with gT\ :sup:`2` being proportional to the wavelength.


.. figure:: ./Mehaute_nonlinear_waves.png
   :width: 100%
   :align: center

This case tests demonstrates the ability of Proteus to simulate the 
generation, propagation and absorption of regular non-linear waves. 

Running the test case
-----

The test case comprises a simple rectangular tank with generation zone at the left side ('x-') and absoprtion zone at the right side ('x+'). To run the test case type:

```
parun regular_waves.py --TwoPhaseFlow -f  -v -D result_folder
```

Wave properties can be modified by the commandline, using for example:

```
parun regular_waves.py --TwoPhaseFlow -f regular_waves.py -v -D result_folder -C "T=2. H=0.05"
```

To run in parallel (example with mpirun and 12 processors):

```
mpirun -np 12 parun regular_waves.py --TwoPhaseFlow -v -D result_folder -C "Tp=2. Hs=0.2"
```


To see guidance on parun options, you can type  

```
parun -h
```

References
----------

- Fenton JD (1988) The numerical solution of steady water wave 
  problems, Comp and Geosc, 14(3), 357-368
  
- Lé Méhauté, B., (1976). “Introduction to Hydrodynamics and water waves”, Springer-Verlag, New York.






