SPARK_DIR=${HOME}/spark-2.0.1
benchmark_dir=${HOME}/pse-hpc/data-analytics-benchmarks/
size=25600000
iter=20
num_cores=288

${SPARK_DIR}/sbin/stop-all.sh
${SPARK_DIR}/sbin/start-all.sh

${SPARK_DIR}/bin/spark-submit $benchmark_dir/src/main/python/logistic_regression_gen.py $size &> tmp_spark_logistic_regression_gen.txt

${SPARK_DIR}/sbin/stop-all.sh

echo -e "-1\n$iter\n$size" > ./logistic_regression_gen.data
MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"

mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./logistic_regression_gen &> tmp_hpat_logistic_regression_gen.txt

