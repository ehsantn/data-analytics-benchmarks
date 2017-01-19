#!/bin/bash -l

#SBATCH --partition=debug
#SBATCH --nodes=64
#SBATCH --time=00:30:00
#SBATCH -e HPAT_res_%j.err
#SBATCH -o HPAT_res_%j.out
#SBATCH -C haswell 

export OMP_NUM_THREADS=1
echo "KMEANS"
srun -n 2048 -c 2 ./kmeans_gen
echo "logistic"
srun -n 2048 -c 2 ./logistic_regression_gen
echo "linear"
srun -n 2048 -c 2 ./linear_regression_gen
echo "kernelscore"
srun -n 2048 -c 2 ./kernelscore_gen

echo "manual KMEANS"
srun -n 2048 -c 2 ./kmeans_manual_gen 320000000 20 5
echo "manual logistic"
srun -n 2048 -c 2 ./logistic_regression_manual_gen 2000000000 20
echo "manual linear"
srun -n 2048 -c 2 ./linear_regression_manual_gen 2000000000 20
echo "manual kernelscore"
srun -n 2048 -c 2 ./kernelscore_gen_manual 2000000000
