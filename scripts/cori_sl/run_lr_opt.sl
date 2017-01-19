#!/bin/bash -l

#SBATCH --partition=debug
#SBATCH --nodes=64
#SBATCH --time=00:30:00
#SBATCH -e HPAT_lropt_%j.err
#SBATCH -o HPAT_lropt_%j.out
#SBATCH -C haswell 

export OMP_NUM_THREADS=1
echo "LR"
srun -n 32 -c 2 ./logistic_regression_gen
srun -n 128 -c 2 ./logistic_regression_gen
srun -n 512 -c 2 ./logistic_regression_gen
srun -n 2048 -c 2 ./logistic_regression_gen

echo "LR noopt"
srun -n 32 -c 2 ./logistic_regression_gen_noopt
srun -n 128 -c 2 ./logistic_regression_gen_noopt
srun -n 512 -c 2 ./logistic_regression_gen_noopt
srun -n 2048 -c 2 ./logistic_regression_gen_noopt

echo "manual LR"
srun -n 32 -c 2 ./logistic_regression_manual_gen 1000000000 20 
srun -n 128 -c 2 ./logistic_regression_manual_gen 1000000000 20
srun -n 512 -c 2 ./logistic_regression_manual_gen 1000000000 20
srun -n 2048 -c 2 ./logistic_regression_manual_gen 1000000000 20
