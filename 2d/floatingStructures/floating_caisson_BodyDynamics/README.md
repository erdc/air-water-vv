# Oscillation of Floating Caisson (2D) [WIP]

## Benchmark Layout

This benchmark consists of testing the roll motion (free or under wave loads) of a floating caisson in 2 dimensions. The computational domain is a rectangular tank with default dimensions of 5m x 1.2m, with default floating caisson dimensions are 0.3m x 0.1m. Initially, water in the tank is at rest and the default water level is 0.9m. The walls of the tank have no slip boundary conditions and the top is left open.
This case works with Body Dynamics.

![Alt text](floating_caisson.png)

## Running
-----

To run the test case type:

```
parun --TwoPhaseFlow -f floating_caisson2D.py -v -D result_folder
```

Geometry and set up options can be modified by the commandline, using for example:

```
parun --TwoPhaseFlow -f floating_caisson2D.py -v -D result_folder -C "T=10."
```

To run in parallel (example with mpirun and 12 processors):

```
mpirun -np 12 parun --TwoPhaseFlow -f floating_caisson2D.py -v -D result_folder -C "T=10."
```


To see guidance on parun options, you can type  

```
parun -h
```


## Context Options


| Options        | Explanation                                                         | Default value |
|----------------|---------------------------------------------------------------------|---------------|
| water_level    | Height of free surface above bottom                                 | 0.9           |
| tank_dim       | Dimensions of the tank                                              | (1., 1.2)     |
| dim            | X-dimension of the caisson2D,Y-dimension of the caisson2D           | (0.3,0.1)     |
| center         | Coordinates of caisson center                                       | (0.5, 0.9)    |
| width          | Width of the caisson                                                | 0.9           |
| free_x         | Translational degrees of freedom                                    | (0., 0., 0.)  |
| free_r         | Rotational degrees of freedom                                       | (0., 0., 1.)  |
| caisson_inertia| Inertia of the caisson [kg.m2]                                      | 0.236         |
