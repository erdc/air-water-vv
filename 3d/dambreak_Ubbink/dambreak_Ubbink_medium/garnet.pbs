#!/bin/bash
#PBS -A ERDCV00898R40
#PBS -l walltime=001:00:00
#PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -q debug
#PBS -N ubbink3dm
#PBS -j oe
#PBS -l application=proteus
#PBS -V
#PBS -m eba
#PBS -M cekees@gmail.com
source /opt/modules/default/etc/modules.sh
source /lustre/shared/projects/proteus/garnet.gnu.sh
cd $PBS_O_WORKDIR
mkdir $WORKDIR/dambreak2D.$PBS_JOBID
aprun -n 32  parun dambreak_Ubbink_medium_so.py -l 5 -v -O ../../../inputTemplates/petsc.options.superlu_dist -D $WORKDIR/dambreak2D.$PBS_JOBID
