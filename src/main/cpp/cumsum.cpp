#include <mpi.h>
#include "hdf5.h"
#include <stdint.h>
#include <float.h>
#include <limits.h>
#include <complex>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include "/home/etotoni/.julia/v0.5/ParallelAccelerator/src/../deps/include/j2c-array.h"
#include "/home/etotoni/.julia/v0.5/ParallelAccelerator/src/../deps/include/pse-types.h"
#include "/home/etotoni/.julia/v0.5/ParallelAccelerator/src/../deps/include/cgen_intrinsics.h"
#include <sstream>
#include <vector>
#include <string>
#include <ctime>
#include <stdlib.h>
#include "/home/etotoni/.julia/v0.5/HPAT/src/../deps/include/hpat.h"
unsigned main_count = 0;
typedef struct
{
} pppcumsum_testp271;
void ppcumsum_testp271(ASCIIString&  file_name, double * __restrict ret0)
{
    pppcumsum_testp271 pselfp;
    j2c_array< int64_t >  _df_id;
    j2c_array< double >  _df_x;
    j2c_array< double >  _df_y;
    j2c_array< double >  cs;
    j2c_array< double >  c;
    int64_t __hpat_h5_dim_size_1_1;
    int64_t SSAValue0;
    hsize_t* SSAValue3;
    double SSAValue1;
    int64_t parfor_index_1_1;
    double parallel_ir_array_temp__6_3_1;
    int64_t parallel_ir_save_array_len_1_1;
    double parallel_ir_reduction_output_1;
    double SSAValue2;
    int32_t __hpat_num_pes;
    int32_t __hpat_node_id;
    int64_t __hpat_dist_arr_start_1;
    int64_t __hpat_dist_arr_div_1;
    int64_t __hpat_dist_arr_count_1;
    int64_t __hpat_dist_arr_start_2;
    int64_t __hpat_dist_arr_div_2;
    int64_t __hpat_dist_arr_count_2;
    int64_t __hpat_dist_arr_start_3;
    int64_t __hpat_dist_arr_div_3;
    int64_t __hpat_dist_arr_count_3;
    int64_t __hpat_dist_arr_start_4;
    int64_t __hpat_dist_arr_div_4;
    int64_t __hpat_dist_arr_count_4;
    int64_t __hpat_loop_start_1;
    int64_t __hpat_loop_end_1;
    int64_t __hpat_loop_div_1;
    double __hpat_reduce_5;
    ;;
    MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
    MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
    hid_t plist_id_1 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_1 != -1);
    herr_t ret_1;
    hid_t file_id_1;
    ret_1 = H5Pset_fapl_mpio(plist_id_1, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_1 != -1);
    file_id_1 = H5Fopen((const char*)file_name.data.data, H5F_ACC_RDONLY, plist_id_1);
    assert(file_id_1 != -1);
    ret_1 = H5Pclose(plist_id_1);
    assert(ret_1 != -1);
    hid_t dataset_id_1;
    dataset_id_1 = H5Dopen2(file_id_1, "/id", H5P_DEFAULT);
    assert(dataset_id_1 != -1);
    ;
    hid_t space_id_1 = H5Dget_space(dataset_id_1);
    assert(space_id_1 != -1);
    hsize_t data_ndim_1 = H5Sget_simple_extent_ndims(space_id_1);
    hsize_t space_dims_1[data_ndim_1];
    H5Sget_simple_extent_dims(space_id_1, space_dims_1, NULL);
    SSAValue3 = space_dims_1;;
    __hpat_h5_dim_size_1_1 = SSAValue3[1-1];
    __hpat_dist_arr_div_1 = (__hpat_h5_dim_size_1_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_1 = (__hpat_node_id) * (__hpat_dist_arr_div_1);
    __hpat_dist_arr_count_1 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_1_1-__hpat_node_id*__hpat_dist_arr_div_1 : __hpat_dist_arr_div_1);
    _df_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
    hsize_t CGen_HDF5_start_1[data_ndim_1];
    hsize_t CGen_HDF5_count_1[data_ndim_1];
    CGen_HDF5_start_1[0] = __hpat_dist_arr_start_1;
    CGen_HDF5_count_1[0] = __hpat_dist_arr_count_1;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_1; i_CGen_dim++)
    {
        CGen_HDF5_start_1[i_CGen_dim] = 0;
        CGen_HDF5_count_1[i_CGen_dim] = space_dims_1[i_CGen_dim];
    }
    ret_1 = H5Sselect_hyperslab(space_id_1, H5S_SELECT_SET, CGen_HDF5_start_1, NULL, CGen_HDF5_count_1, NULL);
    assert(ret_1 != -1);
    hid_t mem_dataspace_1 = H5Screate_simple (data_ndim_1, CGen_HDF5_count_1, NULL);
    assert (mem_dataspace_1 != -1);
    hid_t xfer_plist_1 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_1 != -1);
    double h5_read_start_1 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_1, H5FD_MPIO_COLLECTIVE);
    ret_1 = H5Dread(dataset_id_1, H5T_NATIVE_LLONG, mem_dataspace_1, space_id_1, xfer_plist_1, _df_id.getData());
    assert(ret_1 != -1);
    ;
    ;
    H5Dclose(dataset_id_1);
    H5Fclose(file_id_1);
    ;
    hid_t plist_id_2 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_2 != -1);
    herr_t ret_2;
    hid_t file_id_2;
    ret_2 = H5Pset_fapl_mpio(plist_id_2, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_2 != -1);
    file_id_2 = H5Fopen((const char*)file_name.data.data, H5F_ACC_RDONLY, plist_id_2);
    assert(file_id_2 != -1);
    ret_2 = H5Pclose(plist_id_2);
    assert(ret_2 != -1);
    hid_t dataset_id_2;
    dataset_id_2 = H5Dopen2(file_id_2, "/x", H5P_DEFAULT);
    assert(dataset_id_2 != -1);
    ;
    hid_t space_id_2 = H5Dget_space(dataset_id_2);
    assert(space_id_2 != -1);
    hsize_t data_ndim_2 = H5Sget_simple_extent_ndims(space_id_2);
    hsize_t space_dims_2[data_ndim_2];
    H5Sget_simple_extent_dims(space_id_2, space_dims_2, NULL);
    ;
    __hpat_dist_arr_div_2 = (__hpat_h5_dim_size_1_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_2 = (__hpat_node_id) * (__hpat_dist_arr_div_2);
    __hpat_dist_arr_count_2 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_1_1-__hpat_node_id*__hpat_dist_arr_div_2 : __hpat_dist_arr_div_2);
    _df_x = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_2);
    hsize_t CGen_HDF5_start_2[data_ndim_2];
    hsize_t CGen_HDF5_count_2[data_ndim_2];
    CGen_HDF5_start_2[0] = __hpat_dist_arr_start_2;
    CGen_HDF5_count_2[0] = __hpat_dist_arr_count_2;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_2; i_CGen_dim++)
    {
        CGen_HDF5_start_2[i_CGen_dim] = 0;
        CGen_HDF5_count_2[i_CGen_dim] = space_dims_2[i_CGen_dim];
    }
    ret_2 = H5Sselect_hyperslab(space_id_2, H5S_SELECT_SET, CGen_HDF5_start_2, NULL, CGen_HDF5_count_2, NULL);
    assert(ret_2 != -1);
    hid_t mem_dataspace_2 = H5Screate_simple (data_ndim_2, CGen_HDF5_count_2, NULL);
    assert (mem_dataspace_2 != -1);
    hid_t xfer_plist_2 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_2 != -1);
    double h5_read_start_2 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_2, H5FD_MPIO_COLLECTIVE);
    ret_2 = H5Dread(dataset_id_2, H5T_NATIVE_DOUBLE, mem_dataspace_2, space_id_2, xfer_plist_2, _df_x.getData());
    assert(ret_2 != -1);
    ;
    ;
    H5Dclose(dataset_id_2);
    H5Fclose(file_id_2);
    ;
    hid_t plist_id_3 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_3 != -1);
    herr_t ret_3;
    hid_t file_id_3;
    ret_3 = H5Pset_fapl_mpio(plist_id_3, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_3 != -1);
    file_id_3 = H5Fopen((const char*)file_name.data.data, H5F_ACC_RDONLY, plist_id_3);
    assert(file_id_3 != -1);
    ret_3 = H5Pclose(plist_id_3);
    assert(ret_3 != -1);
    hid_t dataset_id_3;
    dataset_id_3 = H5Dopen2(file_id_3, "/y", H5P_DEFAULT);
    assert(dataset_id_3 != -1);
    ;
    hid_t space_id_3 = H5Dget_space(dataset_id_3);
    assert(space_id_3 != -1);
    hsize_t data_ndim_3 = H5Sget_simple_extent_ndims(space_id_3);
    hsize_t space_dims_3[data_ndim_3];
    H5Sget_simple_extent_dims(space_id_3, space_dims_3, NULL);
    ;
    __hpat_dist_arr_div_3 = (__hpat_h5_dim_size_1_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_3 = (__hpat_node_id) * (__hpat_dist_arr_div_3);
    __hpat_dist_arr_count_3 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_1_1-__hpat_node_id*__hpat_dist_arr_div_3 : __hpat_dist_arr_div_3);
    _df_y = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_3);
    hsize_t CGen_HDF5_start_3[data_ndim_3];
    hsize_t CGen_HDF5_count_3[data_ndim_3];
    CGen_HDF5_start_3[0] = __hpat_dist_arr_start_3;
    CGen_HDF5_count_3[0] = __hpat_dist_arr_count_3;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_3; i_CGen_dim++)
    {
        CGen_HDF5_start_3[i_CGen_dim] = 0;
        CGen_HDF5_count_3[i_CGen_dim] = space_dims_3[i_CGen_dim];
    }
    ret_3 = H5Sselect_hyperslab(space_id_3, H5S_SELECT_SET, CGen_HDF5_start_3, NULL, CGen_HDF5_count_3, NULL);
    assert(ret_3 != -1);
    hid_t mem_dataspace_3 = H5Screate_simple (data_ndim_3, CGen_HDF5_count_3, NULL);
    assert (mem_dataspace_3 != -1);
    hid_t xfer_plist_3 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_3 != -1);
    double h5_read_start_3 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_3, H5FD_MPIO_COLLECTIVE);
    ret_3 = H5Dread(dataset_id_3, H5T_NATIVE_DOUBLE, mem_dataspace_3, space_id_3, xfer_plist_3, _df_y.getData());
    assert(ret_3 != -1);
    ;
    ;
    H5Dclose(dataset_id_3);
    H5Fclose(file_id_3);
    ;
    t1 = MPI_Wtime();
    SSAValue0 = __hpat_h5_dim_size_1_1;
    __hpat_dist_arr_div_4 = (SSAValue0) / (__hpat_num_pes);
    __hpat_dist_arr_start_4 = (__hpat_node_id) * (__hpat_dist_arr_div_4);
    __hpat_dist_arr_count_4 = ((__hpat_node_id==__hpat_num_pes-1) ? SSAValue0-__hpat_node_id*__hpat_dist_arr_div_4 : __hpat_dist_arr_div_4);
    c = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_4);
    cs = c;
    int64_t size__df_x = _df_x.ARRAYSIZE(1);
    double tmp_sum__df_x=0, prefix__df_x=0;
    for(int i=0; i<size__df_x; i++)
        tmp_sum__df_x += _df_x.data[i];
    MPI_Exscan(&tmp_sum__df_x, &prefix__df_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(int i=0; i<size__df_x; i++)
    {
        prefix__df_x += _df_x.data[i];
        c.data[i] = prefix__df_x;
    }
    ;
    t2 = MPI_Wtime();
    parallel_ir_save_array_len_1_1 = SSAValue0;
    parallel_ir_reduction_output_1 = 0.0;
    __hpat_loop_div_1 = (parallel_ir_save_array_len_1_1) / (__hpat_num_pes);
    __hpat_loop_start_1 = ((__hpat_node_id) * (__hpat_loop_div_1)) + (1);
    __hpat_loop_end_1 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_1 : (__hpat_node_id+1)*__hpat_loop_div_1);
    for ( parfor_index_1_1 = __hpat_loop_start_1; parfor_index_1_1 <= (int64_t)__hpat_loop_end_1; parfor_index_1_1 += 1)
    {
        ;
        parallel_ir_array_temp__6_3_1 = cs.ARRAYELEM(((parfor_index_1_1) - (__hpat_loop_start_1)) + (1));
        SSAValue2 = (parallel_ir_reduction_output_1) + (parallel_ir_array_temp__6_3_1);
        parallel_ir_reduction_output_1 = SSAValue2;
    }
    ;
    __hpat_reduce_5 = 0;
    MPI_Allreduce(&parallel_ir_reduction_output_1, &__hpat_reduce_5, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
    parallel_ir_reduction_output_1 = __hpat_reduce_5;
    SSAValue1 = parallel_ir_reduction_output_1;
    *ret0 = SSAValue1;
    return;

}


void ppcumsum_testp271_unaliased(ASCIIString& __restrict file_name, double * __restrict ret0)
{
    pppcumsum_testp271 pselfp;
    j2c_array< int64_t >  _df_id;
    j2c_array< double >  _df_x;
    j2c_array< double >  _df_y;
    j2c_array< double >  cs;
    j2c_array< double >  c;
    int64_t __hpat_h5_dim_size_1_1;
    int64_t SSAValue0;
    hsize_t* SSAValue3;
    double SSAValue1;
    int64_t parfor_index_1_1;
    double parallel_ir_array_temp__6_3_1;
    int64_t parallel_ir_save_array_len_1_1;
    double parallel_ir_reduction_output_1;
    double SSAValue2;
    int32_t __hpat_num_pes;
    int32_t __hpat_node_id;
    int64_t __hpat_dist_arr_start_1;
    int64_t __hpat_dist_arr_div_1;
    int64_t __hpat_dist_arr_count_1;
    int64_t __hpat_dist_arr_start_2;
    int64_t __hpat_dist_arr_div_2;
    int64_t __hpat_dist_arr_count_2;
    int64_t __hpat_dist_arr_start_3;
    int64_t __hpat_dist_arr_div_3;
    int64_t __hpat_dist_arr_count_3;
    int64_t __hpat_dist_arr_start_4;
    int64_t __hpat_dist_arr_div_4;
    int64_t __hpat_dist_arr_count_4;
    int64_t __hpat_loop_start_1;
    int64_t __hpat_loop_end_1;
    int64_t __hpat_loop_div_1;
    double __hpat_reduce_5;
    ;;
    MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
    MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
    hid_t plist_id_1 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_1 != -1);
    herr_t ret_1;
    hid_t file_id_1;
    ret_1 = H5Pset_fapl_mpio(plist_id_1, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_1 != -1);
    file_id_1 = H5Fopen((const char*)file_name.data.data, H5F_ACC_RDONLY, plist_id_1);
    assert(file_id_1 != -1);
    ret_1 = H5Pclose(plist_id_1);
    assert(ret_1 != -1);
    hid_t dataset_id_1;
    dataset_id_1 = H5Dopen2(file_id_1, "/id", H5P_DEFAULT);
    assert(dataset_id_1 != -1);
    ;
    hid_t space_id_1 = H5Dget_space(dataset_id_1);
    assert(space_id_1 != -1);
    hsize_t data_ndim_1 = H5Sget_simple_extent_ndims(space_id_1);
    hsize_t space_dims_1[data_ndim_1];
    H5Sget_simple_extent_dims(space_id_1, space_dims_1, NULL);
    SSAValue3 = space_dims_1;;
    __hpat_h5_dim_size_1_1 = SSAValue3[1-1];
    __hpat_dist_arr_div_1 = (__hpat_h5_dim_size_1_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_1 = (__hpat_node_id) * (__hpat_dist_arr_div_1);
    __hpat_dist_arr_count_1 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_1_1-__hpat_node_id*__hpat_dist_arr_div_1 : __hpat_dist_arr_div_1);
    _df_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
    hsize_t CGen_HDF5_start_1[data_ndim_1];
    hsize_t CGen_HDF5_count_1[data_ndim_1];
    CGen_HDF5_start_1[0] = __hpat_dist_arr_start_1;
    CGen_HDF5_count_1[0] = __hpat_dist_arr_count_1;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_1; i_CGen_dim++)
    {
        CGen_HDF5_start_1[i_CGen_dim] = 0;
        CGen_HDF5_count_1[i_CGen_dim] = space_dims_1[i_CGen_dim];
    }
    ret_1 = H5Sselect_hyperslab(space_id_1, H5S_SELECT_SET, CGen_HDF5_start_1, NULL, CGen_HDF5_count_1, NULL);
    assert(ret_1 != -1);
    hid_t mem_dataspace_1 = H5Screate_simple (data_ndim_1, CGen_HDF5_count_1, NULL);
    assert (mem_dataspace_1 != -1);
    hid_t xfer_plist_1 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_1 != -1);
    double h5_read_start_1 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_1, H5FD_MPIO_COLLECTIVE);
    ret_1 = H5Dread(dataset_id_1, H5T_NATIVE_LLONG, mem_dataspace_1, space_id_1, xfer_plist_1, _df_id.getData());
    assert(ret_1 != -1);
    ;
    ;
    H5Dclose(dataset_id_1);
    H5Fclose(file_id_1);
    ;
    hid_t plist_id_2 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_2 != -1);
    herr_t ret_2;
    hid_t file_id_2;
    ret_2 = H5Pset_fapl_mpio(plist_id_2, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_2 != -1);
    file_id_2 = H5Fopen((const char*)file_name.data.data, H5F_ACC_RDONLY, plist_id_2);
    assert(file_id_2 != -1);
    ret_2 = H5Pclose(plist_id_2);
    assert(ret_2 != -1);
    hid_t dataset_id_2;
    dataset_id_2 = H5Dopen2(file_id_2, "/x", H5P_DEFAULT);
    assert(dataset_id_2 != -1);
    ;
    hid_t space_id_2 = H5Dget_space(dataset_id_2);
    assert(space_id_2 != -1);
    hsize_t data_ndim_2 = H5Sget_simple_extent_ndims(space_id_2);
    hsize_t space_dims_2[data_ndim_2];
    H5Sget_simple_extent_dims(space_id_2, space_dims_2, NULL);
    ;
    __hpat_dist_arr_div_2 = (__hpat_h5_dim_size_1_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_2 = (__hpat_node_id) * (__hpat_dist_arr_div_2);
    __hpat_dist_arr_count_2 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_1_1-__hpat_node_id*__hpat_dist_arr_div_2 : __hpat_dist_arr_div_2);
    _df_x = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_2);
    hsize_t CGen_HDF5_start_2[data_ndim_2];
    hsize_t CGen_HDF5_count_2[data_ndim_2];
    CGen_HDF5_start_2[0] = __hpat_dist_arr_start_2;
    CGen_HDF5_count_2[0] = __hpat_dist_arr_count_2;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_2; i_CGen_dim++)
    {
        CGen_HDF5_start_2[i_CGen_dim] = 0;
        CGen_HDF5_count_2[i_CGen_dim] = space_dims_2[i_CGen_dim];
    }
    ret_2 = H5Sselect_hyperslab(space_id_2, H5S_SELECT_SET, CGen_HDF5_start_2, NULL, CGen_HDF5_count_2, NULL);
    assert(ret_2 != -1);
    hid_t mem_dataspace_2 = H5Screate_simple (data_ndim_2, CGen_HDF5_count_2, NULL);
    assert (mem_dataspace_2 != -1);
    hid_t xfer_plist_2 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_2 != -1);
    double h5_read_start_2 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_2, H5FD_MPIO_COLLECTIVE);
    ret_2 = H5Dread(dataset_id_2, H5T_NATIVE_DOUBLE, mem_dataspace_2, space_id_2, xfer_plist_2, _df_x.getData());
    assert(ret_2 != -1);
    ;
    ;
    H5Dclose(dataset_id_2);
    H5Fclose(file_id_2);
    ;
    hid_t plist_id_3 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_3 != -1);
    herr_t ret_3;
    hid_t file_id_3;
    ret_3 = H5Pset_fapl_mpio(plist_id_3, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_3 != -1);
    file_id_3 = H5Fopen((const char*)file_name.data.data, H5F_ACC_RDONLY, plist_id_3);
    assert(file_id_3 != -1);
    ret_3 = H5Pclose(plist_id_3);
    assert(ret_3 != -1);
    hid_t dataset_id_3;
    dataset_id_3 = H5Dopen2(file_id_3, "/y", H5P_DEFAULT);
    assert(dataset_id_3 != -1);
    ;
    hid_t space_id_3 = H5Dget_space(dataset_id_3);
    assert(space_id_3 != -1);
    hsize_t data_ndim_3 = H5Sget_simple_extent_ndims(space_id_3);
    hsize_t space_dims_3[data_ndim_3];
    H5Sget_simple_extent_dims(space_id_3, space_dims_3, NULL);
    ;
    __hpat_dist_arr_div_3 = (__hpat_h5_dim_size_1_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_3 = (__hpat_node_id) * (__hpat_dist_arr_div_3);
    __hpat_dist_arr_count_3 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_1_1-__hpat_node_id*__hpat_dist_arr_div_3 : __hpat_dist_arr_div_3);
    _df_y = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_3);
    hsize_t CGen_HDF5_start_3[data_ndim_3];
    hsize_t CGen_HDF5_count_3[data_ndim_3];
    CGen_HDF5_start_3[0] = __hpat_dist_arr_start_3;
    CGen_HDF5_count_3[0] = __hpat_dist_arr_count_3;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_3; i_CGen_dim++)
    {
        CGen_HDF5_start_3[i_CGen_dim] = 0;
        CGen_HDF5_count_3[i_CGen_dim] = space_dims_3[i_CGen_dim];
    }
    ret_3 = H5Sselect_hyperslab(space_id_3, H5S_SELECT_SET, CGen_HDF5_start_3, NULL, CGen_HDF5_count_3, NULL);
    assert(ret_3 != -1);
    hid_t mem_dataspace_3 = H5Screate_simple (data_ndim_3, CGen_HDF5_count_3, NULL);
    assert (mem_dataspace_3 != -1);
    hid_t xfer_plist_3 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_3 != -1);
    double h5_read_start_3 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_3, H5FD_MPIO_COLLECTIVE);
    ret_3 = H5Dread(dataset_id_3, H5T_NATIVE_DOUBLE, mem_dataspace_3, space_id_3, xfer_plist_3, _df_y.getData());
    assert(ret_3 != -1);
    ;
    ;
    H5Dclose(dataset_id_3);
    H5Fclose(file_id_3);
    ;
    SSAValue0 = __hpat_h5_dim_size_1_1;
    __hpat_dist_arr_div_4 = (SSAValue0) / (__hpat_num_pes);
    __hpat_dist_arr_start_4 = (__hpat_node_id) * (__hpat_dist_arr_div_4);
    __hpat_dist_arr_count_4 = ((__hpat_node_id==__hpat_num_pes-1) ? SSAValue0-__hpat_node_id*__hpat_dist_arr_div_4 : __hpat_dist_arr_div_4);
    c = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_4);
    cs = c;
    int64_t size__df_x = _df_x.ARRAYSIZE(1);
    double tmp_sum__df_x=0, prefix__df_x=0;
    for(int i=0; i<size__df_x; i++)
        tmp_sum__df_x += _df_x.data[i];
    MPI_Exscan(&tmp_sum__df_x, &prefix__df_x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(int i=0; i<size__df_x; i++)
    {
        prefix__df_x += _df_x.data[i];
        c.data[i] = prefix__df_x;
    }
    ;
    parallel_ir_save_array_len_1_1 = SSAValue0;
    parallel_ir_reduction_output_1 = 0.0;
    __hpat_loop_div_1 = (parallel_ir_save_array_len_1_1) / (__hpat_num_pes);
    __hpat_loop_start_1 = ((__hpat_node_id) * (__hpat_loop_div_1)) + (1);
    __hpat_loop_end_1 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_1 : (__hpat_node_id+1)*__hpat_loop_div_1);
    for ( parfor_index_1_1 = __hpat_loop_start_1; parfor_index_1_1 <= (int64_t)__hpat_loop_end_1; parfor_index_1_1 += 1)
    {
        ;
        parallel_ir_array_temp__6_3_1 = cs.ARRAYELEM(((parfor_index_1_1) - (__hpat_loop_start_1)) + (1));
        SSAValue2 = (parallel_ir_reduction_output_1) + (parallel_ir_array_temp__6_3_1);
        parallel_ir_reduction_output_1 = SSAValue2;
    }
    ;
    __hpat_reduce_5 = 0;
    MPI_Allreduce(&parallel_ir_reduction_output_1, &__hpat_reduce_5, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
    parallel_ir_reduction_output_1 = __hpat_reduce_5;
    SSAValue1 = parallel_ir_reduction_output_1;
    *ret0 = SSAValue1;
    return;

}


extern "C" void _ppcumsum_testp271_unaliased_(int run_where, ASCIIString& __restrict file_name , double* __restrict ret0 , bool genMain = true)
{

    if (genMain)
    {
        ++main_count;
        std::stringstream newMain;
        std::stringstream newMainData;
        std::stringstream newMainSh;
        std::stringstream newMainExe;
        newMain << "main" << main_count << ".cc";
        newMainData << "main" << main_count << ".data";
        newMainSh << "main" << main_count << ".sh";
        newMainExe << "main" << main_count;
        std::cout << "Main will be generated in file " << newMain.str() << std::endl;
        std::cout << "Data for main will be in file " << newMainData.str() << std::endl;
        std::cout << "Script to compile is in " << newMainSh.str() << std::endl;
        std::ofstream mainFileData(newMainData.str(), std::ios::out | std::ios::binary);
        mainFileData << run_where << std::endl;
        mainFileData << file_name << std::endl;
        mainFileData.close();
        std::ofstream mainFile(newMain.str());
        mainFile << "#include \"" << __FILE__ << "\"" << std::endl;
        mainFile << "int main(int argc, char *argv[]) {" << std::endl;
        mainFile << "    MPI_Init(&argc, &argv);" << std::endl;
        mainFile << "    std::ifstream mainFileData(\"" << newMainData.str() << "\", std::ios::in | std::ios::binary);" << std::endl;
        mainFile << "    int runwhere;" << std::endl;
        mainFile << "    mainFileData >> runwhere;" << std::endl;
        mainFile << "    double ret0;" << std::endl;
        mainFile << "    ASCIIString file_name;" << std::endl;
        mainFile << "    mainFileData >> file_name;" << std::endl;
        mainFile << "    _ppcumsum_testp271_unaliased_(runwhere, file_name, &ret0, false);" << std::endl;
        mainFile << "    MPI_Finalize();" << std::endl;
        mainFile << "    return 0;" << std::endl;
        mainFile << "}" << std::endl;
        mainFile.close();
        std::ofstream mainFileSh(newMainSh.str());
        mainFileSh << "#!/bin/sh" << std::endl;
        mainFileSh << "mpiicpc -O3 -std=c++11  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << " -lhdf5  -lm " << std::endl;
        mainFileSh.close();
    }

    ppcumsum_testp271_unaliased(file_name, ret0);
}


extern "C" void _ppcumsum_testp271_(int run_where, ASCIIString&  file_name , double* __restrict ret0 , bool genMain = true)
{

    if (genMain)
    {
        ++main_count;
        std::stringstream newMain;
        std::stringstream newMainData;
        std::stringstream newMainSh;
        std::stringstream newMainExe;
        newMain << "main" << main_count << ".cc";
        newMainData << "main" << main_count << ".data";
        newMainSh << "main" << main_count << ".sh";
        newMainExe << "main" << main_count;
        std::cout << "Main will be generated in file " << newMain.str() << std::endl;
        std::cout << "Data for main will be in file " << newMainData.str() << std::endl;
        std::cout << "Script to compile is in " << newMainSh.str() << std::endl;
        std::ofstream mainFileData(newMainData.str(), std::ios::out | std::ios::binary);
        mainFileData << run_where << std::endl;
        mainFileData << file_name << std::endl;
        mainFileData.close();
        std::ofstream mainFile(newMain.str());
        mainFile << "#include \"" << __FILE__ << "\"" << std::endl;
        mainFile << "int main(int argc, char *argv[]) {" << std::endl;
        mainFile << "    MPI_Init(&argc, &argv);" << std::endl;
        mainFile << "    std::ifstream mainFileData(\"" << newMainData.str() << "\", std::ios::in | std::ios::binary);" << std::endl;
        mainFile << "    int runwhere;" << std::endl;
        mainFile << "    mainFileData >> runwhere;" << std::endl;
        mainFile << "    double ret0;" << std::endl;
        mainFile << "    ASCIIString file_name;" << std::endl;
        mainFile << "    mainFileData >> file_name;" << std::endl;
        mainFile << "    _ppcumsum_testp271_(runwhere, file_name, &ret0, false);" << std::endl;
        mainFile << "    MPI_Finalize();" << std::endl;
        mainFile << "    return 0;" << std::endl;
        mainFile << "}" << std::endl;
        mainFile.close();
        std::ofstream mainFileSh(newMainSh.str());
        mainFileSh << "#!/bin/sh" << std::endl;
        mainFileSh << "mpiicpc -O3 -std=c++11  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << " -lhdf5  -lm " << std::endl;
        mainFileSh.close();
    }

    ppcumsum_testp271(file_name, ret0);
}


extern "C"
void *j2c_array_new(int key, void*data, unsigned ndim, int64_t *dims)
{
    void *a = NULL;
    switch(key)
    {
        case 1:
            a = new  j2c_array< uint8_t > ((uint8_t*)data, ndim, dims);
            break;
        default:
            fprintf(stderr, "j2c_array_new called with invalid key %d", key);
            assert(false);
            break;
    }
    return a;
}
