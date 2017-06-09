# Oscillation of Floating Caisson (3D) [WIP]

## Benchmark Layout

This benchmark consists of testing the roll, pitch or heave motion of a floating caisson. The computational domain is a cuboidal tank with default dimensions of 6m x 1.5m x 1.5m, and the default floating caisson dimensions are 0.594/ x 0.416m x 0.640m. Water in the tank is at rest for initial conditions and the default water level is 0.6m. The walls of the tank have free slip boundary conditions and the top is left open.

## Running

The benchmark can be run using the following command:
```
parun caisson3D_oscillation_so.py -l 2 -v -O petsc_options -D output_folder -C context_options
```
where:
* `petsc_options` must point to the petsc options file
* `output_folder` is the name of the folder for the output files
* `context_options` are options for running the benchmark (see section below)

## Context Options


### Predefined Test Cases

Test cases can be run with the following option:

```
-C "test_case"=4
```
The type of test being run is as follows:
* test_type=1: heave
* test_type=2: pitch
* test_type=3: roll

The possible test cases set the following values automatically:

| Test   | Length | Width | Height |   VCG | Draft |    Ixx |    Iyy |    Izz |
|--------|--------|-------|--------|-------|-------|--------|--------|--------|
| **1**  |  0.594 | 0.416 |  0.640 | 0.191 | 0.425 |  4.956 |  6.620 |  5.515 |
| **2**  |  0.594 | 0.416 |  0.640 | 0.201 | 0.500 |  5.356 |  7.220 |  6.215 |
| **3**  |  0.594 | 0.608 |  0.640 | 0.188 | 0.350 |  9.089 |  8.848 |  9.307 |
| **4**  |  0.594 | 0.608 |  0.640 | 0.189 |       |        |        |        |
| **5**  |  0.594 | 0.608 |  0.640 | 0.189 | 0.425 |  9.789 |  9,448 | 10.807 |
| **6**  |  0.594 | 0.608 |  0.640 | 0.203 | 0.500 | 11.089 | 10.648 | 12.307 |
| **7**  |  1.143 | 0.800 |  0.640 |       |       |        |        |        |
| **8**  |  1.143 | 0.800 |  0.640 | 0.178 | 0.350 | 28.944 | 47.389 | 58.050 |
| **9**  |  1.143 | 0.800 |  0.640 | 0.183 | 0.425 | 32.544 | 54.689 | 68.550 |
| **10** |  1.143 | 0.800 |  0.640 | 0.200 | 0.500 | 37.044 | 62.889 | 79.050 |




### All Options
The default values are equivalent to context options of `"test_case"=1` and `"test_type"=3` 

| Options        | Explanation                                         | Default value       |
|----------------|-----------------------------------------------------|---------------------|
| test_case      | Test case number                                    | None                |
| test_type      | Type of test                                        | 3                   |
| water_level    | Height of free surface above bottom                 | 0.6                 |
| tank_dim       | Dimensions of the tank                              | (5,1.5,1.5)         |
| tank_sponge    | Length of absorption zones (front/back, left/right) | (2.,0.)             |
| caisson_dim    | Dimensions of the floating caisson                  | (0.594,0.416,0.640) |
| free_x         | Translational degrees of freedom                    | (0,0,0)             |
| free_r         | Rotational degrees of freedom                       | (1,0,0)             |
| VCG            | Vertical height of the center of gravity            | 0.175               |
| draft          | Draft of the caisson                                | 0.425               |
| It             | Inertia tensor: Ixx, Iyy, and Izz components        | (4.956,6.620,5.515) |
| rotation_angle | Angle of initial rotation (in radians)              | pi/12               |
| rotation_axis  | Axis of initial rotation                            | (1,0,0)             |
