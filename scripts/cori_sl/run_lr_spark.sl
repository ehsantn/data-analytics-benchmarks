#!/bin/bash

#SBATCH -p debug
#SBATCH -N 16
#SBATCH -t 00:30:00
#SBATCH -e spark_lr_16_5_%j.err
#SBATCH -o spark_lr_16_5_%j.out
#SBATCH -C haswell 

module load spark

which start-all.sh

start-all.sh

spark-submit --conf spark.kryoserializer.buffer.max=1024m ../src/main/python/logistic_regression_gen.py 1000000000 5 10

stop-all.sh

