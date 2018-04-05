#!/bin/bash
#SBATCH -J wavetank	     #Job Name	
#SBATCH -o oe.wavetank.asm.o%j   #Output and Error File
#SBATCH -n 192		     #Number of mpi asks
#SBATCH -p development	     #Queue
#SBATCH -t 01:00:00	     #Run Time
#SBATCH --mail-user=steve.a.mattis@gmail.com	#email
#SBATCH --mail-type=begin				#when to email
#SBATCH --mail-type=end                               #when to email 
#SBATCH -A ADCIRC			#account
set -x
#source ${WORK}/src/proteus-hashstack8/envConfig
#source /scratch/01082/smattis/src/proteus-hashdist-11-2015/envConfig
source /home1/01082/smattis/src/proteus/envConfig
mkdir $SLURM_JOB_NAME.$SLURM_JOB_ID
#/usr/local/bin/ibrun /scratch/01082/smattis/src/proteus-hashdist-11-2015/proteus/stampede.gnu/bin/python /scratch/01082/smattis/src/proteus-hashdist-11-2015/proteus/stampede.gnu/bin/parun tank_so.py  -l 3 -v -D $SLURM_JOB_NAME.$SLURM_JOB_ID -O ../../../inputTemplates/petsc.options.asm -o context.options -p
ibrun parun tank_so.py  -l 3 -v -D $SLURM_JOB_NAME.$SLURM_JOB_ID -O ../../petsc.options.asm -o context.options #-p