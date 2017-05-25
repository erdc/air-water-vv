#!/bin/sh

#PBS -S /bin/sh
#PBS -N cas1215_Sharp_crested_weir_VM_V2
#PBS -o sharpweir.output.o   # stdout file
#PBS -e sharpweir.output.e   # stderr file
#PBS -l nodes=2:ppn=12 # nodes required / processors per node
#PBS -q highp          # queue name

source /etc/profile.d/modules.sh


module load proteus/0.9.0/
cd $PBS_O_WORKDIR

export NCPUS="24"
mpirun -n $NCPUS  parun sharp_crested_weir_so.py -m -p -l 7 -v -O /apps/proteus/0.9.0/proteus-mprans/benchmarks/inputTemplates/petsc.options.superlu_dist -D sharp_crested_weir.$PBS_JOBID | tee log.$PBS_JOBID
