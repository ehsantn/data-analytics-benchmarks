#!/bin/bash -l

#SBATCH --partition=debug
#SBATCH --nodes=64
#SBATCH --time=00:10:00
#SBATCH -e HPAT_res_%j.err
#SBATCH -o HPAT_res_%j.out
#SBATCH -C haswell 

export OMP_NUM_THREADS=1
echo "KMEANS"
srun -n 32 -c 2 ./kmeans_gen
srun -n 128 -c 2 ./kmeans_gen
srun -n 512 -c 2 ./kmeans_gen
srun -n 2048 -c 2 ./kmeans_gen

echo "KMEANS noopt"
srun -n 32 -c 2 ./kmeans_gen_noopt
srun -n 128 -c 2 ./kmeans_gen_noopt
srun -n 512 -c 2 ./kmeans_gen_noopt
srun -n 2048 -c 2 ./kmeans_gen_noopt

echo "manual KMEANS"
srun -n 32 -c 2 ./kmeans_manual_gen 320000000 20 5
srun -n 128 -c 2 ./kmeans_manual_gen 320000000 20 5
srun -n 512 -c 2 ./kmeans_manual_gen 320000000 20 5
srun -n 2048 -c 2 ./kmeans_manual_gen 320000000 20 5
