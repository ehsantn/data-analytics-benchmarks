SPARK_DIR=${HOME}/spark-2.0.1
data_path=${HOME}/tmp/
benchmark_dir=${HOME}/pse-hpc/spark-sql-query-tests/

${SPARK_DIR}/sbin/stop-all.sh
${SPARK_DIR}/sbin/start-all.sh

${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 $benchmark_dir/src/main/python/udf_v1.py $data_path/udf_large.csv &> tmp_spark_udf.txt
${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 $benchmark_dir/src/main/python/udf_v2.py $data_path/udf_large.csv &> tmp_spark_udf2.txt

${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 --class Udf_v1 ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $data_path/udf_large.csv &> tmp_spark_udf3.txt
${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 --class Udf_v2 ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $data_path/udf_large.csv &> tmp_spark_udf4.txt

${SPARK_DIR}/sbin/stop-all.sh

MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"

mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./udf_v1 &> tmp_hpat_udf1.txt
mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./udf_v2 &> tmp_hpat_udf2.txt
