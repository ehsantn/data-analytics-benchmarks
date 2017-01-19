SPARK_DIR=${HOME}/spark-2.0.2
benchmark_dir=${HOME}/data-analytics-benchmarks/
size=256000000
num_cores=288

${SPARK_DIR}/sbin/stop-all.sh
${SPARK_DIR}/sbin/start-all.sh

${SPARK_DIR}/bin/spark-submit --class KernelScore3 ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar $size &> tmp_spark_kernelscore_gen.txt

${SPARK_DIR}/sbin/stop-all.sh

scala -J-Xmx28g -classpath ${benchmark_dir}/target/scala-2.11/benchmarks_2.11-0.3.jar KernelScoreS $size $num_cores &> tmp_kernelscore_scala.txt

julia $benchmark_dir/src/main/julia/kernelscore_gen_jl.jl &> tmp_kernelscore_julia.txt

echo -e "-1\n$size" > ./kernelscore_gen.data
#MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"
#mpiicpc -std=c++11 -O3 -xHost -o kernelscore_gen kernelscore_gen_main.cc
#mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./kernelscore_gen &> tmp_hpat_kernelscore_gen.txt
mpirun -n 32 -ppn 32 ./kernelscore_gen &> tmp_hpat_kernelscore_gen.txt
