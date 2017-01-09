#include <random>
#include <mpi.h>
#include <stdint.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <stdlib.h>
#include <assert.h>
double t1,t2;
double* logistic_regression(int iterations, int D, int64_t N)
{
    double *points, *labels, *w;
    int64_t i;
    int64_t ind;

    int64_t kk;
    int32_t num_pes;
    int32_t node_id;

    int64_t my_count;
    int64_t my_start;
    int64_t pe_div;

    //std::random_device rand_device;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    //std::default_random_engine rand_generator(rand_device());
    std::default_random_engine rand_generator(0);
    MPI_Comm_size(MPI_COMM_WORLD,&num_pes);
    MPI_Comm_rank(MPI_COMM_WORLD,&node_id);
    pe_div = (N) / (num_pes);
    my_start = (node_id) * (pe_div);
    my_count = ((node_id==num_pes-1) ? N-node_id*pe_div : pe_div);
    points = new double[D*my_count];
    labels = new double[my_count];
    for(int64_t i=0; i<my_count*D; i++)
        points[i] = distribution(rand_generator);
    for(int64_t i=0; i<my_count; i++)
        labels[i] = distribution(rand_generator);
    w = new double[D];
    if (node_id==0) {
        for ( i = 0; i < (int64_t)D; i += 1)
        {
            w[i] = 2.0*distribution(rand_generator)-1.0;
        }
    }
    MPI_Bcast(w, D, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* glob_grad = new double[D];
    double* grad = new double[D];
    
    t1 = MPI_Wtime();
    for(int64_t iter=0; iter<iterations; iter++){
        for(int kk=0; kk<D; kk++)
            grad[kk] = 0.0;
        for ( ind = 0; ind < my_count; ind +=1)
        {
                double wpp = 0.0;
                for(int kk=0; kk<D; kk++)
                    wpp += w[kk]*points[ind*D+kk];
                double calc = ((1.0/(1.0+expf(-labels[ind]*wpp)))-1.0)*labels[ind];
                for(int kk=0; kk<D; kk++)
                    grad[kk] += calc*points[ind*D+kk];

        }
        MPI_Allreduce(grad, glob_grad, D, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int kk = 0; kk < D; kk += 1)
        {
            w[kk] -= glob_grad[kk];
        }
    }
    t2 = MPI_Wtime();
    delete[] grad;
    delete[] glob_grad;
    delete[] labels;
    delete[] points;
    return w;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int iter = 5;
    int D = 10;
    int64_t N = 10;
    if(argc>2) {
        N = atoi(argv[1]);
        iter = atoi(argv[2]);
    }
    printf("N %d D %d iter :%d\n", N, D, iter);
    double* res = logistic_regression(iter,D,N);
    int node_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&node_id);
    if (node_id==0) printf("iter:%d N:%ld\n res:%lf %lf exec time:%lf\n", iter,N,res[0],res[1],t2-t1);
    MPI_Finalize();
    return 0;
}

