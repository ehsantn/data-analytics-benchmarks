SPARK_DIR=${HOME}/spark-2.0.2
benchmark_dir=${HOME}/data-analytics-benchmarks/

size=64000000
iter=20
D=10
num_center=5

num_cores=36
total_cores=144 
#total_cores=36

${SPARK_DIR}/sbin/stop-all.sh
${SPARK_DIR}/sbin/start-all.sh

${SPARK_DIR}/bin/spark-submit $benchmark_dir/src/main/python/kmeans_gen.py $size $iter $D $num_center &> res_spark_kmeans_gen.out

${SPARK_DIR}/sbin/stop-all.sh

#HOSTS="ip-172-31-43-11.us-west-2.compute.internal,ip-172-31-43-10.us-west-2.compute.internal,ip-172-31-43-9.us-west-2.compute.internal,ip-172-31-43-8.us-west-2.compute.internal"
#HOSTS=localhost
#mpirun -n 144 -ppn 36 -hosts ${HOSTS} python $benchmark_dir/src/main/python/logistic_regression_mpi.py $size $iter &> tmp_python_kmeans_gen.txt

echo -e "-1\n$num_center\n$iter\n$size" > ./kmeans_gen.data

#MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"
#mpiicpc -std=c++11 -O3 -xHost -o kmeans_gen kmeans_gen_main.cc -mkl
#mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./kmeans_gen &> tmp_hpat_kmeans_gen.txt

mpirun -n $total_cores -ppn $num_cores -hosts ${HOSTS} ./kmeans_gen &> res_hpat_kmeans_gen.out

mpirun -n $total_cores -ppn $num_cores -hosts ${HOSTS} ./kmeans_manual_gen $size $iter $num_center &> res_manual_kmeans_gen.out
