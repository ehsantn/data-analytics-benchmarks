SPARK_DIR=${HOME}/spark-2.0.1
data_path=${HOME}/tmp/
benchmark_dir=${HOME}/pse-hpc/spark-sql-query-tests/
file_name=cumsum_large

${SPARK_DIR}/sbin/stop-all.sh
${SPARK_DIR}/sbin/start-all.sh

${SPARK_DIR}/bin/spark-submit --conf spark.sql.autoBroadcastJoinThreshold=-1 --class SimpleMovingAvg ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $data_path/${file_name}.csv &> tmp_spark_smovingavg.txt

${SPARK_DIR}/sbin/stop-all.sh

python3 $benchmark_dir/src/main/python/smovingavg_pd.py $data_path/${file_name}.csv &> pandas_smovingavg.txt

julia $benchmark_dir/src/main/julia/smovingavg_jldf.jl $data_path/${file_name}.hdf5 &> julia_smovingavg.txt

MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"

mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./smovingavg &> tmp_hpat_smovingavg.txt
