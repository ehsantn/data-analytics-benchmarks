SPARK_DIR=${HOME}/spark-2.0.2
benchmark_dir=${HOME}/data-analytics-benchmarks/

size=256000000
iter=20
D=10
F=4

num_cores=36
#total_cores=144 
total_cores=36

${SPARK_DIR}/sbin/stop-all.sh
${SPARK_DIR}/sbin/start-all.sh

${SPARK_DIR}/bin/spark-submit $benchmark_dir/src/main/python/linear_regression_gen.py $size $iter $D $F &> res_spark_linear_regression_gen.out

${SPARK_DIR}/sbin/stop-all.sh

#HOSTS="ip-172-31-43-11.us-west-2.compute.internal,ip-172-31-43-10.us-west-2.compute.internal,ip-172-31-43-9.us-west-2.compute.internal,ip-172-31-43-8.us-west-2.compute.internal"
#HOSTS=localhost
#mpirun -n 144 -ppn 36 -hosts ${HOSTS} python $benchmark_dir/src/main/python/linear_regression_mpi.py $size $iter &> tmp_python_linear_regression_gen.txt

echo -e "-1\n$iter\n$size" > ./linear_regression_gen.data

#MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"
#mpiicpc -std=c++11 -O3 -xHost -o linear_regression_gen linear_regression_gen_main.cc -mkl
#mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./linear_regression_gen &> tmp_hpat_linear_regression_gen.txt

mpirun -n $total_cores -ppn $num_cores -hosts ${HOSTS} ./linear_regression_gen &> res_hpat_linear_regression_gen.out

mpirun -n $total_cores -ppn $num_cores -hosts ${HOSTS} ./linear_regression_manual_gen $size $iter &> res_manual_linear_regression_gen.out
