#!/bin/bash -l

#SBATCH --partition=debug
#SBATCH --nodes=64
#SBATCH --time=00:30:00
#SBATCH -e HPAT_lropt_%j.err
#SBATCH -o HPAT_lropt_%j.out
#SBATCH -C haswell 

export OMP_NUM_THREADS=1
echo "LR"
srun -N 1 -n 32 -c 2 --cpu_bind=cores ./logistic_regression_gen
srun -N 4 -n 128 -c 2 --cpu_bind=cores ./logistic_regression_gen
srun -N 16 -n 512 -c 2 --cpu_bind=cores ./logistic_regression_gen
srun -N 64 -n 2048 -c 2 --cpu_bind=cores ./logistic_regression_gen

echo "LR noopt"
srun -N 1 -n 32 -c 2 --cpu_bind=cores ./logistic_regression_gen_noopt
srun -N 4 -n 128 -c 2 --cpu_bind=cores ./logistic_regression_gen_noopt
srun -N 16 -n 512 -c 2 --cpu_bind=cores ./logistic_regression_gen_noopt
srun -N 64 -n 2048 -c 2 --cpu_bind=cores ./logistic_regression_gen_noopt

echo "manual LR"
srun -N 1 -n 32 -c 2 --cpu_bind=cores ./logistic_regression_manual_gen 1000000000 20 
srun -N 4 -n 128 -c 2 --cpu_bind=cores ./logistic_regression_manual_gen 1000000000 20
srun -N 16 -n 512 -c 2 --cpu_bind=cores ./logistic_regression_manual_gen 1000000000 20
srun -N 64 -n 2048 -c 2 --cpu_bind=cores ./logistic_regression_manual_gen 1000000000 20
