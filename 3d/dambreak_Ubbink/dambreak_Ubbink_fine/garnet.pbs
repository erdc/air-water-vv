#!/bin/bash
#PBS -A ERDCV00898R40
#PBS -l select=8:ncpus=32:mpiprocs=32
#PBS -l walltime=008:00:00
#PBS -q standard
#PBS -N ubbink3df
#PBS -j oe
#PBS -l application=proteus
#PBS -V
#PBS -m eba
#PBS -M cekees@gmail.com
source /opt/modules/default/etc/modules.sh
source /lustre/shared/projects/proteus/garnet.gnu.sh
cd $PBS_O_WORKDIR
mkdir $WORKDIR/ubbink3df.$PBS_JOBID
aprun -n 256 parun dambreak_Ubbink_fine_so.py -F -l 5 -v -O ../../../inputTemplates/petsc.options.asm -D $WORKDIR/ubbink3df.$PBS_JOBID
