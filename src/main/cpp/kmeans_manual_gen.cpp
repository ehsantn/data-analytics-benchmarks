
#include <random>
#include <mpi.h>
#include "hdf5.h"
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

double* kmeans(int numCenter, int64_t iterations, int D, int64_t N)
{
    double *points, *centroids;
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
    for(int64_t i=0; i<my_count*D; i++)
        points[i] = distribution(rand_generator);
  
    double* centroids_sum = new double[D*numCenter];
    uint64_t* centroids_count = new uint64_t[numCenter];

    double *glob_centroids_sum = new double[D*numCenter];
    uint64_t* glob_centroids_count = new uint64_t[numCenter];

    centroids = new double[D*numCenter];
    if (node_id==0) {
        for ( i = 0; i < D*numCenter; i += 1)
        {
            centroids[i] = distribution(rand_generator);
        }
    }
    MPI_Bcast(centroids, D*numCenter, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for(int64_t iter=0; iter<iterations; iter++){
        // reset centroid sums and counts
        for(int cc=0; cc<numCenter; cc++)
        {
            centroids_count[cc] = 0;
            for(int kk=0; kk<D; kk++)
                centroids_sum[cc*D+kk] = 0.0;
        }

        for ( ind = 0; ind < my_count; ind +=1)
        {
            int best_ind = 0;
            double best_distance = DBL_MAX;
            for(int cc=0; cc<numCenter; cc++)
            {
                double tmpDist = 0.0;
                for(int kk=0; kk<D; kk++)
                {
                    tmpDist += pow(points[ind*D+kk]-centroids[cc*D+kk],2.0);
                }
                tmpDist = sqrt(tmpDist);
                if(tmpDist<best_distance)
                {
                    best_ind = cc;
                    best_distance = tmpDist;
                }
            }
            centroids_count[best_ind]++;
            for(int kk=0; kk<D; kk++)
                centroids_sum[best_ind*D+kk] += points[ind*D+kk];

        }
        MPI_Allreduce(centroids_sum, glob_centroids_sum, D*numCenter, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(centroids_count, glob_centroids_count, numCenter, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        for(int cc=0; cc<numCenter; cc++)
        {
            for (int kk = 0; kk < D; kk += 1)
            {
                centroids[cc*D+kk] = glob_centroids_sum[cc*D+kk]/glob_centroids_count[cc];
            }
        }
    }
    delete[] glob_centroids_sum;
    delete[] glob_centroids_count;
    delete[] points;
    return centroids;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int iter = 5;
    int num_center=5;
    int D = 10;
    int64_t N = 16;
    if(argc>1) {
        N = atoi(argv[1]);
        iter = atoi(argv[2]);
        num_center = atoi(argv[3]);
    }

    double t1 = MPI_Wtime();
    double* res = kmeans(num_center,iter, D, N);
    double t2 = MPI_Wtime();
    int node_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&node_id);
    if (node_id==0) printf("iter:%d N:%lld\n res:%lf %lf exec time:%lf\n", iter,N,res[0],res[1],t2-t1);
    MPI_Finalize();
    return 0;
}

