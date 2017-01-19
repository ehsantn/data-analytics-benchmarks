#!/bin/bash

#SBATCH -p debug
#SBATCH -N 64
#SBATCH -t 00:30:00
#SBATCH -e spark_res_ks_%j.err
#SBATCH -o spark_res_ks_%j.out
#SBATCH -C haswell 

module load spark

which start-all.sh

start-all.sh

spark-submit --conf spark.kryoserializer.buffer.max=1024m ../src/main/python/logistic_regression_gen.py 2000000000 20 10
spark-submit --conf spark.kryoserializer.buffer.max=1024m ../src/main/python/kmeans_gen.py 320000000 20 10 5
spark-submit --conf spark.kryoserializer.buffer.max=1024m ../src/main/python/linear_regression_gen.py 2000000000 20 10 4
spark-submit --conf spark.kryoserializer.buffer.max=1024m ../src/main/python/kernelscore_gen.py 2000000000

stop-all.sh

