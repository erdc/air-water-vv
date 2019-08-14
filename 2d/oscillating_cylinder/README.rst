OSCILLATING CYLINDER 
===================

Description
-----------
::

Predicting oscillatory motion in submerged pipelines is a key application in fluid/structure interaction and a
challenging problem to address with computacional and physical modelling due to the complexity of the processes
involved in coupling highly turbulent flows with the pipeline motion. The experimental data used are whose found 
in Fu et al., (2014)

::
 
According with this experimental configuration, a cylinder with diameter D=0.25 was placed in a 196 m long, 10m wide 
and 4.2m deep towing tank. In the laboratory tests a 2.5m long, 2.4m wide and 0.003m thick steel plate was placed
near the bottom of the flume. The plate could be adjusted to different levels in order to mimic different gap 
heights from the seabed. The cylinder and the plate were towed with a speed of 0.8m/s and a mechanical vertical
oscillatory motion was forced. The scope of the experiment was to explore the interaction of vortex shedding in
different vertical motion periods that could represent i.e. wave-induced motion.

::
 
Proteus is used to simulate this experimental configuration, by forcing a vertical oscillation in the pipeline in the form:


-y(t)= Yo sin(2πfot) 

Where:

-Yo: the amplitude of the oscillation
-fo: the oscillating frequency

At the same time a steady current is applied in the [x-] boundary with a velocity of magnitude 0.8 m/s in the x-axes in a 2D numerical tank.

The geometry of the experiment is illustrated in the following figure.

 
.. figure:: ./oscillating_cylinder_drawing.png



Context Options
---------------

+---------------------+--------------------------------------------------------------+--------------------+
| Options             | Explanation                                                  | Default value      |
+=====================+==============================================================+====================+
| mwl                 | Height of free surface above bottom                          | 2.0                |
+---------------------+--------------------------------------------------------------+--------------------+
| tank_dim            | Dimensions of the tank                                       | (4.75,2.5)         |
+---------------------+--------------------------------------------------------------+--------------------+
| current             | Enabling generation of steady current                        | True               |
+---------------------+--------------------------------------------------------------+--------------------+
| U                   | Steady velocity of the current                               | [0.8,0.,0.]        |
+---------------------+--------------------------------------------------------------+--------------------+
| rampTime            | Duration in which the velocity of the current is established | 10.                |
+---------------------+--------------------------------------------------------------+--------------------+
| circle2D            | Switch on/off the extistance of the pipeline in the domain   | True               |
+---------------------+--------------------------------------------------------------+--------------------+
| circleBC            | Boundary Conditions in the pipe                              | NoSlip             |
+---------------------+--------------------------------------------------------------+--------------------+
| InputMotion         | Forcing Oscillation in the pipe                              | True               |
+---------------------+--------------------------------------------------------------+--------------------+
| At                  | Amplitude of imposed sinusoidal translational motion         | [0.0, 0.075, 0.0]  |
+---------------------+--------------------------------------------------------------+--------------------+
| Tt                  | Period of imposed sinusoidal translational motion            | [0.0, 1.30, 0.0]   |
+---------------------+--------------------------------------------------------------+--------------------+
| refinement_level    | Used to define element size he=radius/refinement_level       | 2.5                |
+---------------------+--------------------------------------------------------------+--------------------+



Running
-----

To run the test case type:

```
parun --TwoPhaseFlow -f oscillating_cylinder.py -v -D result_folder
```

Geometry and setup options can be modified by the commandline, using for example:

```
parun --TwoPhaseFlow -f oscillating_cylinder.py -v -D result_folder -C "U=[0.5,0.,0.]"
```

To run in parallel (example with mpirun and 12 processors):

```
mpirun -np 12 parun --TwoPhaseFlow -f oscillating_cylinder.py -v -D result_folder -C "U=[0.5,0.,0.]"
```


To see guidance on parun options, you can type  

```
parun -h
```

 
References 
----------
* Fu S, Xu Y; and Chen Y (2014), Seabed Effects on the Hydrodinamics of a Circular Cylinder Undergoing 
  Vortex-Induced Vibration at High Reynolds Number, Journal of Waterway, Port, Coastal and Ocean 
  Engineering-ASCE, 140, 04014008
