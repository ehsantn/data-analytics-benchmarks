set -x
SPARK_DIR=${HOME}/spark-2.0.1
data_path=${HOME}/tmp/
benchmark_dir=${HOME}/pse-hpc/spark-sql-query-tests/
file_name=udf_large

${SPARK_DIR}/sbin/stop-all.sh
${SPARK_DIR}/sbin/start-all.sh

${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 $benchmark_dir/src/main/python/aggregate.py $data_path/${file_name}.csv &> tmp_spark_aggregate_p.txt

${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 --class AggregateDF1 ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $data_path/${file_name}.csv &> tmp_spark_aggregate1.txt
${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 --class AggregateDF2 ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $data_path/${file_name}.csv &> tmp_spark_aggregate2.txt

${SPARK_DIR}/sbin/stop-all.sh

python3 $benchmark_dir/src/main/python/aggregate_pd.py $data_path/${file_name}.csv &> pandas_aggregate.txt

julia $benchmark_dir/src/main/julia/aggregate_jldf.jl $data_path/${file_name}.hdf5 &> julia_aggregate.txt


MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"

mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./aggregate &> tmp_hpat_aggregate.txt

