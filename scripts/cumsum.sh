SPARK_DIR=${HOME}/spark-2.0.1
data_path=${HOME}/tmp/
benchmark_dir=${HOME}/pse-hpc/spark-sql-query-tests/

${SPARK_DIR}/sbin/stop-all.sh
${SPARK_DIR}/sbin/start-all.sh

${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 --class Cumsum ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $data_path/cumsum_large.csv &> tmp_spark_cumsum.txt

${SPARK_DIR}/sbin/stop-all.sh

MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"

mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./cumsum &> tmp_hpat_cumsum.txt
