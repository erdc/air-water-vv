#!/bin/bash

parun marin_so.py -l 5 -v -C "T=2.0 schur_solver='two_phase_PCD' ns_forceStrongDirichlet=True stabilization='proteus_full' A_block_amg=True Refinement=10 useOnlyVF=True" -O ../../inputTemplates/petsc.options.schur.pcd.amg -D tests/pcd_r5_mpi2_A_amg_lu
