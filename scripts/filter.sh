SPARK_DIR=${HOME}/spark-2.0.1
data_path=${HOME}/tmp/
benchmark_dir=${HOME}/pse-hpc/spark-sql-query-tests/
file_name=filter_2b

${SPARK_DIR}/sbin/stop-all.sh
${SPARK_DIR}/sbin/start-all.sh

${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 $benchmark_dir/src/main/python/filter.py $data_path/${file_name}.csv &> tmp_spark_filter_p.txt

${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 --class FilterDF1 ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $data_path/${file_name}.csv &> tmp_spark_filter1.txt
${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 --class FilterDF2 ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $data_path/${file_name}.csv &> tmp_spark_filter2.txt

${SPARK_DIR}/sbin/stop-all.sh

python3 $benchmark_dir/src/main/python/filter_pd.py $data_path/${file_name}.csv &> pandas_filter.txt

julia $benchmark_dir/src/main/julia/filter_jldf.jl $data_path/${file_name}.hdf5 &> julia_filter.txt


MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"

mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./filter &> tmp_hpat_filter.txt

