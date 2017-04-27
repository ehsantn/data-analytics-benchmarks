#!/bin/bash -l

#SBATCH --partition=debug
#SBATCH --nodes=64
#SBATCH --time=00:30:00
#SBATCH -e py_lr_%j.err
#SBATCH -o py_lr_%j.out
#SBATCH -C haswell 

export OMP_NUM_THREADS=1
srun -N 1 -n 32 -c 2 --cpu_bind=cores python /global/homes/t/totoni/data-analytics-benchmarks/src/main/python/logistic_regression_mpi.py 1000000000 20 
srun -N 4 -n 128 -c 2 --cpu_bind=cores python /global/homes/t/totoni/data-analytics-benchmarks/src/main/python/logistic_regression_mpi.py 1000000000 20
srun -N 16 -n 512 -c 2 --cpu_bind=cores python /global/homes/t/totoni/data-analytics-benchmarks/src/main/python/logistic_regression_mpi.py 1000000000 20
srun -N 64 -n 2048 -c 2 --cpu_bind=cores python /global/homes/t/totoni/data-analytics-benchmarks/src/main/python/logistic_regression_mpi.py 1000000000 20
