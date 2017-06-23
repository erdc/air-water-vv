# Oscillation of Floating Caisson (2D) [WIP]

## Benchmark Layout

This benchmark consists of testing the roll motion (free or under wave loads) of a floating caisson in 2 dimensions. The computational domain is a rectangular tank with default dimensions of 5m x 1.2m, with default floating caisson dimensions are 0.3m x 0.1m. Initially, water in the tank is at rest and the default water level is 0.9m. The walls of the tank have no slip boundary conditions and the top is left open.

![Alt text](floating_caisson.png)

## Running

The benchmark can be run using the following command:
```
parun caisson2D_oscillation_so.py -l 2 -v -O petsc_options -D output_folder -C context_options
```
where:
* `petsc_options` must point to the petsc options file
* `output_folder` is the name of the folder for the output files
* `context_options` are options for running the benchmark (see section below)

## Context Options


| Options        | Explanation                                                         | Default value |
|----------------|---------------------------------------------------------------------|---------------|
| water_level    | Height of free surface above bottom                                 | 0.9           |
| tank_dim       | Dimensions of the tank                                              | (1., 1.2)     |
| tank_sponge    | Length of absorption zones (left, right)                            | (2., 2.)      |
| dimx           | X-dimension of the caisson2D                                        | 0.3           |
| dimy           | Y-dimension of the caisson2D                                        | 0.1           |
| center         | Coordinates of caisson center                                       | (0.5, 0.9)    |
| width          | Width of the caisson                                                | 0.9           |
| free_x         | Translational degrees of freedom                                    | (0., 0., 0.)  |
| free_r         | Rotational degrees of freedom                                       | (0., 0., 1.)  |
| caisson_inertia| Inertia of the caisson [kg.m2]                                      | 0.236         |
| rotation_angle | Angle of initial rotation (in degrees)                              | 15.           |
