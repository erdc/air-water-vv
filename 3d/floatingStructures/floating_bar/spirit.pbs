#!/bin/bash
#PBS -A ERDCV00898R40
##PBS -l walltime=037:00:00
##PBS -l ncpus=1536
#PBS -l walltime=008:00:00
#PBS -l ncpus=16
#PBS -q standard
#PBS -N bar
#PBS -j oe
#PBS -l application=proteus
#PBS -V
#PBS -m eba
#PBS -M cekees@gmail.com
cd $PBS_O_WORKDIR
cp *.py $JOBDIR
cp $HOME/air-water-vv/inputTemplates/petsc.options.asm $JOBDIR
cp $HOME/air-water-vv/inputTemplates/petsc.options.superlu_dist $JOBDIR
cd $JOBDIR
mpiexec_mpt -np ${BC_MPI_TASKS_ALLOC} parun floating_bar_so.py -l 7 -F -C "gen_mesh=True refinement_level=1 parallel=True bar_height=0.6 cfl=0.33 nsave=250 dt_init=0.001 free_r=(1.0,1.0,0.0)" -O petsc.options.asm

