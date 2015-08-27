#!/bin/bash
#SBATCH --comment "GMRES/ls/Arnoldi"
#SBATCH -J " GMRES alone "

#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#SBATCH --time=00:10:00

#SBATCH -c 1
#SBATCH -N 1

#SBATCH --gres=gpu:1

echo $CUDA_VISIBLE_DEVICES
# Retourne 0 ou 1


srun ./a.out --matrix-from-file rdb968.mtx --tolerance 1e-10 --restart 500
