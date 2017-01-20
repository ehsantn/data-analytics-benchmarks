SPARK_DIR=${HOME}/spark-2.0.2
benchmark_dir=${HOME}/data-analytics-benchmarks/

size=256000000
iter=20

num_cores=36
total_cores=144 
#total_cores=36

#${SPARK_DIR}/sbin/stop-all.sh
#${SPARK_DIR}/sbin/start-all.sh

#${SPARK_DIR}/bin/spark-submit --master spark://ip-172-31-20-24.us-west-2.compute.internal:7080 $benchmark_dir/src/main/python/logistic_regression_gen.py $size $iter &> res_spark_logistic_regression_gen.out

#${SPARK_DIR}/sbin/stop-all.sh

HOSTS="ip-172-31-20-24.us-west-2.compute.internal,ip-172-31-18-133.us-west-2.compute.internal,ip-172-31-23-50.us-west-2.compute.internal,ip-172-31-23-37.us-west-2.compute.internal"
#HOSTS=localhost
#mpirun -n 144 -ppn 36 -hosts ${HOSTS} python $benchmark_dir/src/main/python/logistic_regression_mpi.py $size $iter &> tmp_python_logistic_regression_gen.txt

echo -e "-1\n$iter\n$size" > ./logistic_regression_gen.data

#MPI_CONF="-genv I_MPI_ADJUST_ALLGATHERV 4  -genv I_MPI_ADJUST_ALLTOALLV 1 -genv I_MPI_FABRICS shm:dapl"
#mpiicpc -std=c++11 -O3 -xHost -o logistic_regression_gen logistic_regression_gen_main.cc -mkl
#mpirun ${MPI_CONF} -hosts psephi07-ib,psephi08-ib,psephi09-ib,psephi10-ib -n 144 -ppn 36 ./logistic_regression_gen &> tmp_hpat_logistic_regression_gen.txt

mpirun -n $total_cores -ppn $num_cores -hosts ${HOSTS} ./logistic_regression_gen &> res_hpat_logistic_regression_gen.out

mpirun -n $total_cores -ppn $num_cores -hosts ${HOSTS} ./logistic_regression_manual_gen $size $iter &> res_manual_logistic_regression_gen.out
