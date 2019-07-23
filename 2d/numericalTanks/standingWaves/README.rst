Standing Waves
====================================

Description
----------

Standing waves are formed when regular waves interact with a fully reflective wall. In this case, a distinctive patterns emerges and the wave amplitude becomes a sinusioidal function in space. The reflected waves are fully synconised with the incident ones at an even number of half wavelenghts from the wall, whilst cancelling each other at an odd number of wavelength, this creating the distinctive wave patterns of nodes and antinodes, where the amplitude is double and zero, respectively.

This test case is used to demonstrate the capability of Proteus of modelling standing wave patterns and in particular of absorbing the reflected waves at the generation / absorption zone.

Test case
-----

The test case comprises a simple rectangular tank with generation zone at the left side ('x-') and a reflective wall at the right side ('x+'). To run the test case type:

```
parun --TwoPhaseFlow -f standing_waves.py -v -D result_folder
```

Wave properties can be modified by the commandline, using for example:

```
parun --TwoPhaseFlow -f standing_waves.py -v -C -D result_folder "T=2 H=0.05"
```

To run in parallel (example with mpirun and 12 processors):

```
mpirun -np 12 parun --TwoPhaseFlow -f standing_waves.py -v -C -D result_folder "Tp=2 Hs=0.2"
```


To see guidance on parun options, you can type  

```
parun -h
```







