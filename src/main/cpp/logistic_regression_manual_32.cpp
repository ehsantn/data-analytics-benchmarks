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

float* logistic_regression(int64_t iterations, char* file_name)
{
    float *points, *labels, *w;
    int64_t i;
    int64_t ind;
    int64_t N;

    int64_t kk;
    int64_t D;
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
    hid_t plist_id_1 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_1 != -1);
    herr_t ret_1;
    hid_t file_id_1;
    ret_1 = H5Pset_fapl_mpio(plist_id_1, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_1 != -1);
    file_id_1 = H5Fopen((const char*)file_name, H5F_ACC_RDONLY, plist_id_1);
    assert(file_id_1 != -1);
    ret_1 = H5Pclose(plist_id_1);
    assert(ret_1 != -1);
    hid_t dataset_id_1;
    dataset_id_1 = H5Dopen2(file_id_1, "/points", H5P_DEFAULT);
    assert(dataset_id_1 != -1);
    hid_t space_id_1 = H5Dget_space(dataset_id_1);
    assert(space_id_1 != -1);
    hsize_t data_ndim_1 = H5Sget_simple_extent_ndims(space_id_1);
    hsize_t space_dims_1[data_ndim_1];
    H5Sget_simple_extent_dims(space_id_1, space_dims_1, NULL);
    D = space_dims_1[1];
    N = space_dims_1[0]; 
    pe_div = (N) / (num_pes);
    my_start = (node_id) * (pe_div);
    my_count = ((node_id==num_pes-1) ? N-node_id*pe_div : pe_div);
    points = new float[D*my_count];
    hsize_t HDF5_start_1[data_ndim_1];
    hsize_t HDF5_count_1[data_ndim_1];
    HDF5_start_1[0] = my_start;
    HDF5_count_1[0] = my_count;
    HDF5_start_1[1] = 0;
    HDF5_count_1[1] = space_dims_1[1];
    ret_1 = H5Sselect_hyperslab(space_id_1, H5S_SELECT_SET, HDF5_start_1, NULL, HDF5_count_1, NULL);
    assert(ret_1 != -1);
    hid_t mem_dataspace_1 = H5Screate_simple (data_ndim_1, HDF5_count_1, NULL);
    assert (mem_dataspace_1 != -1);
    hid_t xfer_plist_1 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_1 != -1);
    H5Pset_dxpl_mpio(xfer_plist_1, H5FD_MPIO_COLLECTIVE);
    ret_1 = H5Dread(dataset_id_1, H5T_NATIVE_FLOAT, mem_dataspace_1, space_id_1, xfer_plist_1, points);
    assert(ret_1 != -1);
    H5Dclose(dataset_id_1);
    H5Fclose(file_id_1);
    hid_t plist_id_2 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_2 != -1);
    herr_t ret_2;
    hid_t file_id_2;
    ret_2 = H5Pset_fapl_mpio(plist_id_2, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_2 != -1);
    file_id_2 = H5Fopen((const char*)file_name, H5F_ACC_RDONLY, plist_id_2);
    assert(file_id_2 != -1);
    ret_2 = H5Pclose(plist_id_2);
    assert(ret_2 != -1);
    hid_t dataset_id_2;
    dataset_id_2 = H5Dopen2(file_id_2, "/responses", H5P_DEFAULT);
    assert(dataset_id_2 != -1);
    hid_t space_id_2 = H5Dget_space(dataset_id_2);
    assert(space_id_2 != -1);
    hsize_t data_ndim_2 = H5Sget_simple_extent_ndims(space_id_2);
    hsize_t space_dims_2[data_ndim_2];
    H5Sget_simple_extent_dims(space_id_2, space_dims_2, NULL);
    assert(N == space_dims_2[0]);
    labels = new float[my_count];
    hsize_t HDF5_start_2[data_ndim_2];
    hsize_t HDF5_count_2[data_ndim_2];
    HDF5_start_2[0] = my_start;
    HDF5_count_2[0] = my_count;
    ret_2 = H5Sselect_hyperslab(space_id_2, H5S_SELECT_SET, HDF5_start_2, NULL, HDF5_count_2, NULL);
    assert(ret_2 != -1);
    hid_t mem_dataspace_2 = H5Screate_simple (data_ndim_2, HDF5_count_2, NULL);
    assert (mem_dataspace_2 != -1);
    hid_t xfer_plist_2 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_2 != -1);
    H5Pset_dxpl_mpio(xfer_plist_2, H5FD_MPIO_COLLECTIVE);
    ret_2 = H5Dread(dataset_id_2, H5T_NATIVE_FLOAT, mem_dataspace_2, space_id_2, xfer_plist_2, labels);
    assert(ret_2 != -1);
    H5Dclose(dataset_id_2);
    H5Fclose(file_id_2);
    
    w = new float[D];
    if (node_id==0) {
        for ( i = 0; i < (int64_t)D; i += 1)
        {
            w[i] = 2.0*distribution(rand_generator)-1.0;
        }
    }
    MPI_Bcast(w, D, MPI_FLOAT, 0, MPI_COMM_WORLD);

    float* glob_grad = new float[D];
    float* grad = new float[D];
    
    for(int64_t iter=0; iter<iterations; iter++){
        for(int kk=0; kk<D; kk++)
            grad[kk] = 0.0;
        for ( ind = 0; ind < my_count; ind +=1)
        {
                float wpp = 0.0;
                for(int kk=0; kk<D; kk++)
                    wpp += w[kk]*points[ind*D+kk];
                float calc = ((1.0/(1.0+expf(-labels[ind]*wpp)))-1.0)*labels[ind];
                for(int kk=0; kk<D; kk++)
                    grad[kk] += calc*points[ind*D+kk];

        }
        MPI_Allreduce(grad, glob_grad, D, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        for (int kk = 0; kk < D; kk += 1)
        {
            w[kk] -= glob_grad[kk];
        }
    }
    delete[] grad;
    delete[] glob_grad;
    delete[] labels;
    delete[] points;
    return w;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int iter = 5;
    char* file_name = "1D_large.hdf5";
    if(argc>1) {
        iter = atoi(argv[1]);
        file_name = argv[2];
    }

    double t1 = MPI_Wtime();
    float* res = logistic_regression(iter, file_name);
    double t2 = MPI_Wtime();
    int node_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&node_id);
    if (node_id==0) printf("iter:%d file:%s\n res:%lf %lf exec time:%lf\n", iter,file_name,res[0],res[1],t2-t1);
    MPI_Finalize();
    return 0;
}

