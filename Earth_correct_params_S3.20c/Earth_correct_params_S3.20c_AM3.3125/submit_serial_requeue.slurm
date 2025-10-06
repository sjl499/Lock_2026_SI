#!/bin/bash

#SBATCH -p serial_requeue  #send to stewart via general partition
#SBATCH -t 5000
#SBATCH -J EHX001E #name of job
#SBATCH -n 1  #number of cores
#SBATCH -o slurm.out
#SBATCH -e slurm.err
#SBATCH -vvv
#SBATCH --mail-user=slock@fas.harvard.edu #email this user
#SBATCH --mem=200 #memory requested for the whole job


./HERCULESv1.0