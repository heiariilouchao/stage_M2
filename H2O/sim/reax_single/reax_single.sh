#!/bin/sh

#SBATCH --job-name="reax_single"
#SBATCH --partition=umformation
#SBATCH -n 1
#SBATCH --output="output.out"
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user="heiarii.lou-chao@etu.umontpellier.fr"

module purge
module load cv-standard
module load python/2.7.12
module load gcc/4.9.3
module load intel/compiler/64/2017.1.132
module load intel/mpi/64/2017.1.132

mpirun -np 1 ~/scratch/lmp -log data/log.setup -in in.reax_single
