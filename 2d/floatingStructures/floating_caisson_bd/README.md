# Oscillation of Floating Caisson (2D) [WIP]

## Benchmark Layout

This benchmark consists of testing the roll motion (free or under wave loads) of a floating caisson in 2 dimensions. The computational domain is a rectangular tank with default dimensions of 5m x 1.2m, with default floating caisson dimensions are 0.3m x 0.1m. Initially, water in the tank is at rest and the default water level is 0.9m. The walls of the tank have no slip boundary conditions and the top is left open.
This case works with Body Dynamics.

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
| water_level    | Height of free surface above bottom                                 | 0.6           |
| tank_dim       | Dimensions of the tank                                              | (5., 1.2)     |
| tank_sponge    | Length of absorption zones (left, right)                            | (2., 2.)      |
| waves          | Boolean to indicate if waves will be generated                      | False         |
| wave_period    | Period of the waves (if any)                                        | 0.8           |
| wave_height    | Height of the waves (if any)                                        | 0.029         |
| wave_dir       | Direction of the waves (if any)                                     | (1., 0., 0.)  |
| caisson_dim    | Dimensions of the floating caisson                                  | (0.3, 0.1)    |
| caisson_coords | Coordinates of caisson (default: middle of the tank at water level) | None          |
| caisson_width  | Width of the caisson                                                | 0.9           |
| free_x         | Translational degrees of freedom                                    | (0, 0, 0)     |
| free_r         | Rotational degrees of freedom                                       | (0, 0, 1)     |
| VCG            | Height of the center of gravity (default: half height of caisson)   | None          |
| inertia        | Inertia of the caisson                                              | 0.236         |
| rotation_angle | Angle of initial rotation (in radians)                              | pi/12         |
