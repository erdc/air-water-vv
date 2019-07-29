Overtopping over constant slope dike 
==============================================

Description
-----------
This application has been set up to calibrate and evaluate the ability of proteus to calculate overtopping over a constant slope, impermeable dike 
The boundary areas [x+] and [y-] are defined as free slip walls, the upper boundary [y+] is left open and the [x-] is the velocity inlet.

The areas of the geometry are the generation and absorption zones, the wave propagation zone before the obstacle, 
the collection tank and a drainage pipe to ensure that the MWL in the landward tank is not decreasing due to the volume of water that overtops. 

The obstacle can be described as a constant-slope, positive freeboard impermeable dike with crest height Rc=0.1m and slope tanθ=1/4. 
The leeward side of the obstacle is designed as a vertical wall, with zero crest width (sharp-crested structure). 
This case study corresponds to the geometry of one of the tests encountered in the CLASH database (EuroTop 2018) so there are experimental data available for comparison. 

.. figure:: ./new_flume.jpg
   :width: 90%
   :align: center


| Tank Options     | Explanation                                                             | Default value |
|------------------|-------------------------------------------------------------------------|---------------|
| tank_height      | Positive Y-Coordinate of the tank                                       | 0.8           |
| L_back           | Horizontal Dimension of the collection tank after the drainage outlet   | 0.8           |
| tank_depth       | Negative Y-Coordinate of the tank / Lower level of drainage pipe        | 0.5           |
| obs_depth        | Negative Y-Coordinate of the tank / Upper level of drainage pipe        | 0.4           |
| tube             | Drainge pipe diameter                                                   | 0.1           |
| deposit width    | Width of the ovetopping collection tank                                 | (0.5, 0.9)    |



Test case
-----

To run the test case type:

```
parun --TwoPhaseFlow -f Overtopping_constant_slope.py -v -D result_folder
```

Wave properties can be modified by the commandline, using for example:

```
parun --TwoPhaseFlow -f Overtopping_constant_slope.py -v -D result_folder -C "mwl=0.3"
```

To run in parallel (example with mpirun and 12 processors):

```
mpirun -np 12 parun --TwoPhaseFlow -f Overtopping_constant_slope.py -v -D result_folder -C "mwl=0.3"
```


To see guidance on parun options, you can type  

```
parun -h
```


References
----------
EurOtop, 2018.  Manual on wave overtopping of sea defences and related structures.  An overtopping manual largely based on European research, but for worldwide application.  Van der Meer, J.W., Allsop, N.W.H., Bruce, T., De Rouck, J., Kortenhaus, A., Pullen, T., Schüttrumpf, H., Troch, P. and Zanuttigh, B., www.overtopping-manual.com

