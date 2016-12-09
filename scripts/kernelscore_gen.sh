SPARK_DIR=${HOME}/spark-2.0.1
benchmark_dir=${HOME}/pse-hpc/data-analytics-benchmarks/
size=256000000

${SPARK_DIR}/sbin/stop-all.sh
${SPARK_DIR}/sbin/start-all.sh

${SPARK_DIR}/bin/spark-submit --class KernelScore3 ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $size &> tmp_spark_kernelscore_gen.txt

${SPARK_DIR}/sbin/stop-all.sh

scala -J-Xmx28g -classpath ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar KernelScoreS $size

echo -e "-1\n$size" > ./kernelscore_gen.data
MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"

mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./kernelscore_gen &> tmp_hpat_kernelscore_gen.txt

