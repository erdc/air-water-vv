Sharp crested weir benchmark case
=================================

A rectangular sharp-crested weir is a hydraulic structure mainly used as a flow control and 
measurement device in open-channels. It is a channel obstruction made from a thin plate 
with a sharp edge at the top (crest), spanned across the width of the channel. This artificial 
barrier obstructs the natural flow of the water, leading to an increase in the water level 
upstream, before spilling over the structure.
The computational domain is a 2D rectangular box with height equal to 1.8m and the length might 
vary from case to case. The weir has a height of P=1.0m and a small width of b=0.01m.
A uniform velocity distribution from bottom to the free water level is imposed in the left wall 
boundary condition. The top of the domain is left open and the right wall allows 
the flow to leave the domain. 
In the following figure a simple sketch of the structure is presented showing the main parameters.

.. figure:: ./SharpWeir.bmp

where, uo is the approach velocity, H is the upstream potential head, hv is the upstream velocity head, Ht = H + hv is the upstream total head, h is the thickness of the nappe, d1 is the backwater depth beneath the nappe, d2 is the backwater depth downstream of the nappe.

This case tests the ability of PROTEUS to simulate the free-surface evolution and the 
flow separation. The results of the simulations can be compared with the data in the following references.
For more details, see runfiles or references.

References
--------------------------------

- Montes, J.S. (1992). "Curvature Analysis of Spillway Profiles." Proc. 11th Australasian Fluid Mechanics Conference AFMC, Hobart, Australia, Paper 7E-7,2, 941-944.
- U. S. Army Engineer Waterways Experiment Station (WES, 1977)  - HYDRAULIC DESIGN CRITERIA - SHEETS 111-11 to 111-14/1 - Overflow Spillway Crest – Upper Nappe Profile.
