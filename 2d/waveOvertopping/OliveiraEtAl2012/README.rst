Overtopping
====================

This case is simulating the generation and propagation of a regular wave group, which interacts with a levee - type structure and overtops (see figure below). This case was studied experimentally and numerically by Oliveira et al 2012. Waves are generated at the left side of the domain by using a numerical moving paddle in order to exactly reproduce the experiment. The time history of paddle displacement is shown below and it is a very good match of the one shown in Oliveira et al. 2012 

The domain is fitted with an absorption zone behind the paddle to prevent energy build-up. The paddle operates for 4 wave periods, with the paddle motion being ramped up and down for one period. This case study serves to validate the ability of numerical models to predict overtopping. In the current set-up, the mesh size is relatively coarse as this serves as a quick test to verify the implementation of the moving paddle module (among others) and to avoid regression. To achieve a good agreement with the experimental data, mesh size has to be set <1cm.

!(./Structure.png)



!(./Paddle.png)
   


References
----------
Tiago C. A. Oliveira, Agustın Sanchez-Arcilla and Xavier Gironella (2012). Simulation of Wave Overtopping of Maritime Structures in a Numerical Wave Flume, J of Appl. Math. 35 2012, Article ID 246146, 19 pages, doi:10.1155/2012/246146 
