SPARK_DIR=${HOME}/spark-2.0.2
benchmark_dir=${HOME}/data-analytics-benchmarks/

size=256000000

num_cores=36
total_cores=144 
#total_cores=36

#${SPARK_DIR}/sbin/stop-all.sh
#${SPARK_DIR}/sbin/start-all.sh

#${SPARK_DIR}/bin/spark-submit $benchmark_dir/src/main/python/kernelscore_gen.py $size &> res_spark_kernelscore_gen.out

#${SPARK_DIR}/sbin/stop-all.sh

HOSTS="ip-172-31-20-24.us-west-2.compute.internal,ip-172-31-18-133.us-west-2.compute.internal,ip-172-31-23-50.us-west-2.compute.internal,ip-172-31-23-37.us-west-2.compute.internal"
#HOSTS=localhost
#mpirun -n 144 -ppn 36 -hosts ${HOSTS} python $benchmark_dir/src/main/python/kernelscore_mpi.py $size $iter &> tmp_python_kernelscore_gen.txt

echo -e "-1\n$size" > ./kernelscore_gen.data

#MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"
#mpiicpc -std=c++11 -O3 -xHost -o kernelscore_gen kernelscore_gen_main.cc -mkl
#mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./kernelscore_gen &> tmp_hpat_kernelscore_gen.txt

mpirun -n $total_cores -ppn $num_cores -hosts ${HOSTS} ./kernelscore_gen &> res_hpat_kernelscore_gen.out

mpirun -n $total_cores -ppn $num_cores -hosts ${HOSTS} ./kernelscore_gen_manual $size &> res_manual_kernelscore_gen.out
