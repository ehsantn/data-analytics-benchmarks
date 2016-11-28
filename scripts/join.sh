SPARK_DIR=${HOME}/spark-2.0.1
data_path=${HOME}/tmp/
benchmark_dir=${HOME}/pse-hpc/data-analytics-benchmarks/
file_name=join_500k

${SPARK_DIR}/sbin/stop-all.sh
${SPARK_DIR}/sbin/start-all.sh

${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 $benchmark_dir/src/main/python/join.py $data_path/${file_name}.csv &> tmp_spark_join_p.txt

${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 --class JoinDF1 ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $data_path/${file_name}.csv &> tmp_spark_join1.txt
${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 --class JoinDF2 ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $data_path/${file_name}.csv &> tmp_spark_join2.txt

${SPARK_DIR}/sbin/stop-all.sh

python3 $benchmark_dir/src/main/python/join_pd.py $data_path/${file_name}.csv &> pandas_join.txt

julia $benchmark_dir/src/main/julia/join_jldf.jl $data_path/${file_name}.hdf5 &> julia_join.txt


MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"

mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./join &> tmp_hpat_join.txt

