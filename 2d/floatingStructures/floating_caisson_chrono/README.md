# Oscillation of Floating Caisson (2D) [WIP]

## Benchmark Layout

This benchmark consists of testing the roll motion (free or under wave loads) of a floating caisson in 2 dimensions. The computational domain is a rectangular tank with default dimensions of 5m x 1.2m, with default floating caisson dimensions are 0.3m x 0.1m. Initially, water in the tank is at rest and the default water level is 0.9m. The walls of the tank have no slip boundary conditions and the top is left open. This case is using Chrono.

![Alt text](floating_caisson.png)

## Test

The python test file named test_caissonChono.py is made up of two tests:

* The first one is to know if the case can run.
* The second test is to validate the results comparing them to the theory. For this case we will compare the numerical and theoretical period of the rotation angle.

One can run this test file typing py.test --boxed -v test_caissonChrono.py.

## Running

The benchmark can be run using the following command:
```
parun caisson2D_oscillation_so.py -l 2 -v -O petsc_options -D output_folder -C context_options
```
where:
* `petsc_options` must point to the petsc options file
* `output_folder` is the name of the folder for the output files
* `context_options` are options for running the benchmark (see section below)

The benchmark works only using the [mbd/chrono](https://github.com/erdc-cm/proteus/tree/mbd/chrono) proteus branch.

## Context Options


| Options        | Explanation                                                         | Default value |
|----------------|---------------------------------------------------------------------|---------------|
| water_level    | Height of free surface above bottom                                 | 0.9           |
| tank_dim       | Dimensions of the tank                                              | (1., 1.2)     |
| tank_sponge    | Length of absorption zones (left, right)                            | (2., 2.)      |
| waves          | Boolean to indicate if waves will be generated                      | False         |
| caisson_dim    | Dimensions of the floating caisson                                  | (0.3, 0.1)    |
| caisson_coords | Coordinates of caisson (default: middle of the tank at water level) | None          |
| caisson_width  | Width of the caisson                                                | 0.9           |
| free_x         | Translational degrees of freedom                                    | (0., 0., 0.)  |
| free_r         | Rotational degrees of freedom                                       | (0., 0., 1.)  |
| inertia        | Inertia of the caisson                                              | 0.236         |
| rotation_angle | Angle of initial rotation (in degrees)                              | 15.           |
