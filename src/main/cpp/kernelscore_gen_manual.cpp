#include <random>
#include <mpi.h>
#include <stdint.h>
#include <float.h>
#include <limits.h>
#include <complex>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>


double kernel_score(double*, int64_t n);
double* gen_data(int64_t n);

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int64_t n = 100;
    if(argc>1)
        n = atoi(argv[1]);
    int node_id, pes;
    MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
    MPI_Comm_size(MPI_COMM_WORLD, &pes);
    int64_t div = n/pes;
    int64_t chunk = (node_id==pes-1) ? n-node_id*div : div;
    double out;
    printf("generating data: %lld\n", chunk);
    double* X = gen_data(chunk);
    double t1 = MPI_Wtime();
    out = kernel_score(X, chunk);
    double t2 = MPI_Wtime();
    if (node_id==0) printf("res:%lf exec time: %lf\n",out,t2-t1);
    MPI_Finalize();
    return 0;
}

double* gen_data(int64_t n) {
    double* X = new double[n];
    std::random_device cgen_rand_device;
    std::uniform_real_distribution<double> cgen_distribution(0.0,1.0);
    std::default_random_engine cgen_rand_generator(cgen_rand_device());
    for(int i=0; i<n; i++)
        X[i] = cgen_distribution(cgen_rand_generator);
    return X;
}

double kernel_score(double* X, int64_t n) {
    double res = 0.0;
    double b22 = 0.5*0.5*2;
    double lbn = log(3*0.5);

    for(int i=0; i<n; i++) {
        double d1 = -pow(X[i]-(-1.0),2)/b22;
        double d2 = -pow(X[i]-(2.0),2)/b22;
        double d3 = -pow(X[i]-(5.0),2)/b22;
        double m = std::min(d1, std::min(d2,d3));
        double s = exp(d1-m)+exp(d2-m)+exp(d3-m);
        res += m-lbn+log(s);
    }

    return res;
}

