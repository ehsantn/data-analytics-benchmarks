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

#include "/home/etotoni/.julia/v0.5/HPAT/src/../deps/include/hpat_sort.h"
unsigned main_count = 0;
typedef struct
{
} pppcumsum_testp271;
void ppcumsum_testp271(ASCIIString&  file_name, ASCIIString&  file_name2, double * __restrict ret0)
{
    pppcumsum_testp271 pselfp;
    j2c_array< int64_t >  _df1_id;
    j2c_array< double >  _df1_x;
    j2c_array< double >  _df1_y;
    j2c_array< int64_t >  _df2_id;
    j2c_array< double >  _df2_x1;
    j2c_array< double >  _df2_y1;
    j2c_array< int64_t >  _df3_id;
    j2c_array< double >  _df3_x;
    j2c_array< double >  _df3_y;
    j2c_array< double >  _df3_x1;
    j2c_array< double >  _df3_y1;
    int64_t __hpat_h5_dim_size_2_1;
    int64_t __hpat_h5_dim_size_5_1;
    hsize_t* SSAValue0;
    hsize_t* SSAValue1;
    double SSAValue2;
    int64_t parfor_index_1_1;
    double parallel_ir_array_temp__20_3_1;
    int64_t parallel_ir_save_array_len_1_1;
    double parallel_ir_reduction_output_1;
    double SSAValue3;
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
    int64_t __hpat_dist_arr_start_5;
    int64_t __hpat_dist_arr_div_5;
    int64_t __hpat_dist_arr_count_5;
    int64_t __hpat_dist_arr_start_6;
    int64_t __hpat_dist_arr_div_6;
    int64_t __hpat_dist_arr_count_6;
    ;;
    MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
    MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
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
    dataset_id_2 = H5Dopen2(file_id_2, "/id", H5P_DEFAULT);
    assert(dataset_id_2 != -1);
    ;
    hid_t space_id_2 = H5Dget_space(dataset_id_2);
    assert(space_id_2 != -1);
    hsize_t data_ndim_2 = H5Sget_simple_extent_ndims(space_id_2);
    hsize_t space_dims_2[data_ndim_2];
    H5Sget_simple_extent_dims(space_id_2, space_dims_2, NULL);
    SSAValue0 = space_dims_2;;
    __hpat_h5_dim_size_2_1 = SSAValue0[1-1];
    __hpat_dist_arr_div_1 = (__hpat_h5_dim_size_2_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_1 = (__hpat_node_id) * (__hpat_dist_arr_div_1);
    __hpat_dist_arr_count_1 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_2_1-__hpat_node_id*__hpat_dist_arr_div_1 : __hpat_dist_arr_div_1);
    _df1_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
    hsize_t CGen_HDF5_start_2[data_ndim_2];
    hsize_t CGen_HDF5_count_2[data_ndim_2];
    CGen_HDF5_start_2[0] = __hpat_dist_arr_start_1;
    CGen_HDF5_count_2[0] = __hpat_dist_arr_count_1;
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
    ret_2 = H5Dread(dataset_id_2, H5T_NATIVE_LLONG, mem_dataspace_2, space_id_2, xfer_plist_2, _df1_id.getData());
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
    dataset_id_3 = H5Dopen2(file_id_3, "/x", H5P_DEFAULT);
    assert(dataset_id_3 != -1);
    ;
    hid_t space_id_3 = H5Dget_space(dataset_id_3);
    assert(space_id_3 != -1);
    hsize_t data_ndim_3 = H5Sget_simple_extent_ndims(space_id_3);
    hsize_t space_dims_3[data_ndim_3];
    H5Sget_simple_extent_dims(space_id_3, space_dims_3, NULL);
    ;
    __hpat_dist_arr_div_2 = (__hpat_h5_dim_size_2_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_2 = (__hpat_node_id) * (__hpat_dist_arr_div_2);
    __hpat_dist_arr_count_2 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_2_1-__hpat_node_id*__hpat_dist_arr_div_2 : __hpat_dist_arr_div_2);
    _df1_x = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_2);
    hsize_t CGen_HDF5_start_3[data_ndim_3];
    hsize_t CGen_HDF5_count_3[data_ndim_3];
    CGen_HDF5_start_3[0] = __hpat_dist_arr_start_2;
    CGen_HDF5_count_3[0] = __hpat_dist_arr_count_2;
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
    ret_3 = H5Dread(dataset_id_3, H5T_NATIVE_DOUBLE, mem_dataspace_3, space_id_3, xfer_plist_3, _df1_x.getData());
    assert(ret_3 != -1);
    ;
    ;
    H5Dclose(dataset_id_3);
    H5Fclose(file_id_3);
    ;
    hid_t plist_id_4 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_4 != -1);
    herr_t ret_4;
    hid_t file_id_4;
    ret_4 = H5Pset_fapl_mpio(plist_id_4, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_4 != -1);
    file_id_4 = H5Fopen((const char*)file_name.data.data, H5F_ACC_RDONLY, plist_id_4);
    assert(file_id_4 != -1);
    ret_4 = H5Pclose(plist_id_4);
    assert(ret_4 != -1);
    hid_t dataset_id_4;
    dataset_id_4 = H5Dopen2(file_id_4, "/y", H5P_DEFAULT);
    assert(dataset_id_4 != -1);
    ;
    hid_t space_id_4 = H5Dget_space(dataset_id_4);
    assert(space_id_4 != -1);
    hsize_t data_ndim_4 = H5Sget_simple_extent_ndims(space_id_4);
    hsize_t space_dims_4[data_ndim_4];
    H5Sget_simple_extent_dims(space_id_4, space_dims_4, NULL);
    ;
    __hpat_dist_arr_div_3 = (__hpat_h5_dim_size_2_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_3 = (__hpat_node_id) * (__hpat_dist_arr_div_3);
    __hpat_dist_arr_count_3 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_2_1-__hpat_node_id*__hpat_dist_arr_div_3 : __hpat_dist_arr_div_3);
    _df1_y = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_3);
    hsize_t CGen_HDF5_start_4[data_ndim_4];
    hsize_t CGen_HDF5_count_4[data_ndim_4];
    CGen_HDF5_start_4[0] = __hpat_dist_arr_start_3;
    CGen_HDF5_count_4[0] = __hpat_dist_arr_count_3;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_4; i_CGen_dim++)
    {
        CGen_HDF5_start_4[i_CGen_dim] = 0;
        CGen_HDF5_count_4[i_CGen_dim] = space_dims_4[i_CGen_dim];
    }
    ret_4 = H5Sselect_hyperslab(space_id_4, H5S_SELECT_SET, CGen_HDF5_start_4, NULL, CGen_HDF5_count_4, NULL);
    assert(ret_4 != -1);
    hid_t mem_dataspace_4 = H5Screate_simple (data_ndim_4, CGen_HDF5_count_4, NULL);
    assert (mem_dataspace_4 != -1);
    hid_t xfer_plist_4 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_4 != -1);
    double h5_read_start_4 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_4, H5FD_MPIO_COLLECTIVE);
    ret_4 = H5Dread(dataset_id_4, H5T_NATIVE_DOUBLE, mem_dataspace_4, space_id_4, xfer_plist_4, _df1_y.getData());
    assert(ret_4 != -1);
    ;
    ;
    H5Dclose(dataset_id_4);
    H5Fclose(file_id_4);
    ;
    hid_t plist_id_5 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_5 != -1);
    herr_t ret_5;
    hid_t file_id_5;
    ret_5 = H5Pset_fapl_mpio(plist_id_5, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_5 != -1);
    file_id_5 = H5Fopen((const char*)file_name2.data.data, H5F_ACC_RDONLY, plist_id_5);
    assert(file_id_5 != -1);
    ret_5 = H5Pclose(plist_id_5);
    assert(ret_5 != -1);
    hid_t dataset_id_5;
    dataset_id_5 = H5Dopen2(file_id_5, "/id", H5P_DEFAULT);
    assert(dataset_id_5 != -1);
    ;
    hid_t space_id_5 = H5Dget_space(dataset_id_5);
    assert(space_id_5 != -1);
    hsize_t data_ndim_5 = H5Sget_simple_extent_ndims(space_id_5);
    hsize_t space_dims_5[data_ndim_5];
    H5Sget_simple_extent_dims(space_id_5, space_dims_5, NULL);
    SSAValue1 = space_dims_5;;
    __hpat_h5_dim_size_5_1 = SSAValue1[1-1];
    __hpat_dist_arr_div_4 = (__hpat_h5_dim_size_5_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_4 = (__hpat_node_id) * (__hpat_dist_arr_div_4);
    __hpat_dist_arr_count_4 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_5_1-__hpat_node_id*__hpat_dist_arr_div_4 : __hpat_dist_arr_div_4);
    _df2_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_4);
    hsize_t CGen_HDF5_start_5[data_ndim_5];
    hsize_t CGen_HDF5_count_5[data_ndim_5];
    CGen_HDF5_start_5[0] = __hpat_dist_arr_start_4;
    CGen_HDF5_count_5[0] = __hpat_dist_arr_count_4;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_5; i_CGen_dim++)
    {
        CGen_HDF5_start_5[i_CGen_dim] = 0;
        CGen_HDF5_count_5[i_CGen_dim] = space_dims_5[i_CGen_dim];
    }
    ret_5 = H5Sselect_hyperslab(space_id_5, H5S_SELECT_SET, CGen_HDF5_start_5, NULL, CGen_HDF5_count_5, NULL);
    assert(ret_5 != -1);
    hid_t mem_dataspace_5 = H5Screate_simple (data_ndim_5, CGen_HDF5_count_5, NULL);
    assert (mem_dataspace_5 != -1);
    hid_t xfer_plist_5 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_5 != -1);
    double h5_read_start_5 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_5, H5FD_MPIO_COLLECTIVE);
    ret_5 = H5Dread(dataset_id_5, H5T_NATIVE_LLONG, mem_dataspace_5, space_id_5, xfer_plist_5, _df2_id.getData());
    assert(ret_5 != -1);
    ;
    ;
    H5Dclose(dataset_id_5);
    H5Fclose(file_id_5);
    ;
    hid_t plist_id_6 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_6 != -1);
    herr_t ret_6;
    hid_t file_id_6;
    ret_6 = H5Pset_fapl_mpio(plist_id_6, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_6 != -1);
    file_id_6 = H5Fopen((const char*)file_name2.data.data, H5F_ACC_RDONLY, plist_id_6);
    assert(file_id_6 != -1);
    ret_6 = H5Pclose(plist_id_6);
    assert(ret_6 != -1);
    hid_t dataset_id_6;
    dataset_id_6 = H5Dopen2(file_id_6, "/x1", H5P_DEFAULT);
    assert(dataset_id_6 != -1);
    ;
    hid_t space_id_6 = H5Dget_space(dataset_id_6);
    assert(space_id_6 != -1);
    hsize_t data_ndim_6 = H5Sget_simple_extent_ndims(space_id_6);
    hsize_t space_dims_6[data_ndim_6];
    H5Sget_simple_extent_dims(space_id_6, space_dims_6, NULL);
    ;
    __hpat_dist_arr_div_5 = (__hpat_h5_dim_size_5_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_5 = (__hpat_node_id) * (__hpat_dist_arr_div_5);
    __hpat_dist_arr_count_5 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_5_1-__hpat_node_id*__hpat_dist_arr_div_5 : __hpat_dist_arr_div_5);
    _df2_x1 = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_5);
    hsize_t CGen_HDF5_start_6[data_ndim_6];
    hsize_t CGen_HDF5_count_6[data_ndim_6];
    CGen_HDF5_start_6[0] = __hpat_dist_arr_start_5;
    CGen_HDF5_count_6[0] = __hpat_dist_arr_count_5;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_6; i_CGen_dim++)
    {
        CGen_HDF5_start_6[i_CGen_dim] = 0;
        CGen_HDF5_count_6[i_CGen_dim] = space_dims_6[i_CGen_dim];
    }
    ret_6 = H5Sselect_hyperslab(space_id_6, H5S_SELECT_SET, CGen_HDF5_start_6, NULL, CGen_HDF5_count_6, NULL);
    assert(ret_6 != -1);
    hid_t mem_dataspace_6 = H5Screate_simple (data_ndim_6, CGen_HDF5_count_6, NULL);
    assert (mem_dataspace_6 != -1);
    hid_t xfer_plist_6 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_6 != -1);
    double h5_read_start_6 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_6, H5FD_MPIO_COLLECTIVE);
    ret_6 = H5Dread(dataset_id_6, H5T_NATIVE_DOUBLE, mem_dataspace_6, space_id_6, xfer_plist_6, _df2_x1.getData());
    assert(ret_6 != -1);
    ;
    ;
    H5Dclose(dataset_id_6);
    H5Fclose(file_id_6);
    ;
    hid_t plist_id_7 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_7 != -1);
    herr_t ret_7;
    hid_t file_id_7;
    ret_7 = H5Pset_fapl_mpio(plist_id_7, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_7 != -1);
    file_id_7 = H5Fopen((const char*)file_name2.data.data, H5F_ACC_RDONLY, plist_id_7);
    assert(file_id_7 != -1);
    ret_7 = H5Pclose(plist_id_7);
    assert(ret_7 != -1);
    hid_t dataset_id_7;
    dataset_id_7 = H5Dopen2(file_id_7, "/y1", H5P_DEFAULT);
    assert(dataset_id_7 != -1);
    ;
    hid_t space_id_7 = H5Dget_space(dataset_id_7);
    assert(space_id_7 != -1);
    hsize_t data_ndim_7 = H5Sget_simple_extent_ndims(space_id_7);
    hsize_t space_dims_7[data_ndim_7];
    H5Sget_simple_extent_dims(space_id_7, space_dims_7, NULL);
    ;
    __hpat_dist_arr_div_6 = (__hpat_h5_dim_size_5_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_6 = (__hpat_node_id) * (__hpat_dist_arr_div_6);
    __hpat_dist_arr_count_6 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_5_1-__hpat_node_id*__hpat_dist_arr_div_6 : __hpat_dist_arr_div_6);
    _df2_y1 = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_6);
    hsize_t CGen_HDF5_start_7[data_ndim_7];
    hsize_t CGen_HDF5_count_7[data_ndim_7];
    CGen_HDF5_start_7[0] = __hpat_dist_arr_start_6;
    CGen_HDF5_count_7[0] = __hpat_dist_arr_count_6;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_7; i_CGen_dim++)
    {
        CGen_HDF5_start_7[i_CGen_dim] = 0;
        CGen_HDF5_count_7[i_CGen_dim] = space_dims_7[i_CGen_dim];
    }
    ret_7 = H5Sselect_hyperslab(space_id_7, H5S_SELECT_SET, CGen_HDF5_start_7, NULL, CGen_HDF5_count_7, NULL);
    assert(ret_7 != -1);
    hid_t mem_dataspace_7 = H5Screate_simple (data_ndim_7, CGen_HDF5_count_7, NULL);
    assert (mem_dataspace_7 != -1);
    hid_t xfer_plist_7 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_7 != -1);
    double h5_read_start_7 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_7, H5FD_MPIO_COLLECTIVE);
    ret_7 = H5Dread(dataset_id_7, H5T_NATIVE_DOUBLE, mem_dataspace_7, space_id_7, xfer_plist_7, _df2_y1.getData());
    assert(ret_7 != -1);
    ;
    ;
    H5Dclose(dataset_id_7);
    H5Fclose(file_id_7);
    ;
    t1 = MPI_Wtime();
    int join_num_pes_1;
    MPI_Comm join_comm_1;
    join_num_pes_1 = __hpat_num_pes ;
    join_comm_1 = MPI_COMM_WORLD ;
    int * scount_t1_1;
    int * scount_t2_1;
    int * scount_t1_tmp_1;
    int * scount_t2_tmp_1;
    scount_t1_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    memset (scount_t1_1, 0, sizeof(int)* join_num_pes_1);
    scount_t1_tmp_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    memset (scount_t1_tmp_1, 0, sizeof(int)* join_num_pes_1);
    scount_t2_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    memset (scount_t2_1, 0, sizeof(int)* join_num_pes_1);
    scount_t2_tmp_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    memset (scount_t2_tmp_1, 0, sizeof(int)* join_num_pes_1);
    int rsize_t1_1 = 0;
    int rsize_t2_1 = 0;
    int * rcount_t1_1;
    int * rcount_t2_1;
    rcount_t1_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    rcount_t2_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    int * sdis_t1_1;
    int * rdis_t1_1;
    sdis_t1_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    rdis_t1_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    int *sdis_t2_1;
    int *rdis_t2_1;
    sdis_t2_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    rdis_t2_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    int t1c1_length_join_1 = _df1_id.ARRAYLEN() ;
    int t2c1_length_join_1 = _df2_id.ARRAYLEN() ;
    for (int i = 1 ; i <  t1c1_length_join_1 + 1 ; i++)
    {
        int node_id = _df1_id.ARRAYELEM(i) % join_num_pes_1 ;
        scount_t1_1[node_id]++;
    }
    sdis_t1_1[0]=0;
    for(int i=1;i < join_num_pes_1 ;i++)
    {
        sdis_t1_1[i]=scount_t1_1[i-1] + sdis_t1_1[i-1];
    }
    MPI_Alltoall(scount_t1_1,1,MPI_INT,rcount_t1_1,1,MPI_INT, join_comm_1);
    j2c_array< int64_t > _df1_id_tmp_join_1 = j2c_array<int64_t>::new_j2c_array_1d(NULL, t1c1_length_join_1 );
    j2c_array< double > _df1_x_tmp_join_1 = j2c_array<double>::new_j2c_array_1d(NULL, t1c1_length_join_1 );
    j2c_array< double > _df1_y_tmp_join_1 = j2c_array<double>::new_j2c_array_1d(NULL, t1c1_length_join_1 );
    for (int i = 1 ; i <  t1c1_length_join_1 + 1 ; i++)
    {
        int node_id = _df1_id.ARRAYELEM(i) % join_num_pes_1;
        _df1_id_tmp_join_1.ARRAYELEM(sdis_t1_1[node_id]+scount_t1_tmp_1[node_id]+1) = _df1_id.ARRAYELEM(i);
        _df1_x_tmp_join_1.ARRAYELEM(sdis_t1_1[node_id]+scount_t1_tmp_1[node_id]+1) = _df1_x.ARRAYELEM(i);
        _df1_y_tmp_join_1.ARRAYELEM(sdis_t1_1[node_id]+scount_t1_tmp_1[node_id]+1) = _df1_y.ARRAYELEM(i);
        scount_t1_tmp_1[node_id]++;
    }
    for (int i = 1 ; i <  t2c1_length_join_1 + 1 ; i++)
    {
        int node_id = _df2_id.ARRAYELEM(i) % join_num_pes_1 ;
        scount_t2_1[node_id]++;
    }
    sdis_t2_1[0]=0;
    for(int i=1;i < __hpat_num_pes;i++)
    {
        sdis_t2_1[i]=scount_t2_1[i-1] + sdis_t2_1[i-1];
    }
    MPI_Alltoall(scount_t2_1,1,MPI_INT, rcount_t2_1,1,MPI_INT,MPI_COMM_WORLD);
    j2c_array< int64_t > _df2_id_tmp_join_1 = j2c_array<int64_t>::new_j2c_array_1d(NULL, t2c1_length_join_1);
    j2c_array< double > _df2_x1_tmp_join_1 = j2c_array<double>::new_j2c_array_1d(NULL, t2c1_length_join_1);
    j2c_array< double > _df2_y1_tmp_join_1 = j2c_array<double>::new_j2c_array_1d(NULL, t2c1_length_join_1);
    for (int i = 1 ; i <   t2c1_length_join_1 + 1 ; i++)
    {
        int node_id = _df2_id.ARRAYELEM(i) % join_num_pes_1 ;
        _df2_id_tmp_join_1.ARRAYELEM(sdis_t2_1[node_id]+scount_t2_tmp_1[node_id]+1) = _df2_id.ARRAYELEM(i);
        _df2_x1_tmp_join_1.ARRAYELEM(sdis_t2_1[node_id]+scount_t2_tmp_1[node_id]+1) = _df2_x1.ARRAYELEM(i);
        _df2_y1_tmp_join_1.ARRAYELEM(sdis_t2_1[node_id]+scount_t2_tmp_1[node_id]+1) = _df2_y1.ARRAYELEM(i);
        scount_t2_tmp_1[node_id]++;
    }
    rdis_t1_1[0]=0;
    rdis_t2_1[0]=0;
    for(int i=1;i < join_num_pes_1 ;i++)
    {
        rdis_t1_1[i] = rcount_t1_1[i-1] + rdis_t1_1[i-1];
        rdis_t2_1[i] = rcount_t2_1[i-1] + rdis_t2_1[i-1];
    }
    for(int i=0;i< join_num_pes_1 ;i++)
    {
        rsize_t1_1 = rsize_t1_1 + rcount_t1_1[i];
        rsize_t2_1 = rsize_t2_1 + rcount_t2_1[i];
    }
    j2c_array< int64_t > rbuf__df1_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, rsize_t1_1);
    MPI_Alltoallv(_df1_id_tmp_join_1.getData(), scount_t1_1, sdis_t1_1, MPI_INT64_T,
        rbuf__df1_id.getData(), rcount_t1_1, rdis_t1_1, MPI_INT64_T, join_comm_1);
    _df1_id = rbuf__df1_id;
    j2c_array< double > rbuf__df1_x = j2c_array<double>::new_j2c_array_1d(NULL, rsize_t1_1);
    MPI_Alltoallv(_df1_x_tmp_join_1.getData(), scount_t1_1, sdis_t1_1, MPI_INT64_T,
        rbuf__df1_x.getData(), rcount_t1_1, rdis_t1_1, MPI_INT64_T, join_comm_1);
    _df1_x = rbuf__df1_x;
    j2c_array< double > rbuf__df1_y = j2c_array<double>::new_j2c_array_1d(NULL, rsize_t1_1);
    MPI_Alltoallv(_df1_y_tmp_join_1.getData(), scount_t1_1, sdis_t1_1, MPI_INT64_T,
        rbuf__df1_y.getData(), rcount_t1_1, rdis_t1_1, MPI_INT64_T, join_comm_1);
    _df1_y = rbuf__df1_y;
    j2c_array< int64_t > rbuf__df2_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, rsize_t2_1);
    MPI_Alltoallv(_df2_id_tmp_join_1.getData(), scount_t2_1, sdis_t2_1, MPI_INT64_T,
        rbuf__df2_id.getData(), rcount_t2_1, rdis_t2_1, MPI_INT64_T, join_comm_1);
    _df2_id = rbuf__df2_id;
    j2c_array< double > rbuf__df2_x1 = j2c_array<double>::new_j2c_array_1d(NULL, rsize_t2_1);
    MPI_Alltoallv(_df2_x1_tmp_join_1.getData(), scount_t2_1, sdis_t2_1, MPI_INT64_T,
        rbuf__df2_x1.getData(), rcount_t2_1, rdis_t2_1, MPI_INT64_T, join_comm_1);
    _df2_x1 = rbuf__df2_x1;
    j2c_array< double > rbuf__df2_y1 = j2c_array<double>::new_j2c_array_1d(NULL, rsize_t2_1);
    MPI_Alltoallv(_df2_y1_tmp_join_1.getData(), scount_t2_1, sdis_t2_1, MPI_INT64_T,
        rbuf__df2_y1.getData(), rcount_t2_1, rdis_t2_1, MPI_INT64_T, join_comm_1);
    _df2_y1 = rbuf__df2_y1;
    int table_new_counter_join_1 = 1 ;
    std::vector< int64_t> *vec_1__df3_id = new std::vector< int64_t>();
    std::vector< double> *vec_1__df3_x = new std::vector< double>();
    std::vector< double> *vec_1__df3_y = new std::vector< double>();
    std::vector< double> *vec_1__df3_x1 = new std::vector< double>();
    std::vector< double> *vec_1__df3_y1 = new std::vector< double>();
    int64_t * t1_all_arrays_1[3 - 1];
    int64_t * t2_all_arrays_1[3 - 1];
    t1_all_arrays_1 [0] = ( int64_t *) _df1_x.getData();
    t1_all_arrays_1 [1] = ( int64_t *) _df1_y.getData();
    t2_all_arrays_1 [0] = ( int64_t *) _df2_x1.getData();
    t2_all_arrays_1 [1] = ( int64_t *) _df2_y1.getData();
    __hpat_timsort(( int64_t *) _df1_id.getData(), rsize_t1_1 , t1_all_arrays_1, 3 - 1);
    __hpat_timsort(( int64_t *) _df2_id.getData(), rsize_t2_1 , t2_all_arrays_1, 3 - 1);
    int left_join_table_1 = 1;
    int right_join_table_1 = 1;
    while ( (left_join_table_1 < rsize_t1_1 + 1) && (right_join_table_1 < rsize_t2_1 + 1) )
    {
        if(_df1_id.ARRAYELEM(left_join_table_1) == _df2_id.ARRAYELEM(right_join_table_1))
        {
            vec_1__df3_id->push_back( _df1_id.ARRAYELEM(left_join_table_1) );
            vec_1__df3_x->push_back( _df1_x.ARRAYELEM(left_join_table_1) );
            vec_1__df3_y->push_back( _df1_y.ARRAYELEM(left_join_table_1) );
            vec_1__df3_x1->push_back( _df2_x1.ARRAYELEM(right_join_table_1) );
            vec_1__df3_y1->push_back( _df2_y1.ARRAYELEM(right_join_table_1) );
            table_new_counter_join_1++;
            int tmp_left_join_table_1 = left_join_table_1 + 1 ;
            while((tmp_left_join_table_1 < rsize_t1_1 + 1) && (_df1_id.ARRAYELEM(tmp_left_join_table_1) == _df2_id.ARRAYELEM(right_join_table_1)))
            {
                vec_1__df3_id->push_back( _df1_id.ARRAYELEM(tmp_left_join_table_1) );
                vec_1__df3_x->push_back( _df1_x.ARRAYELEM(tmp_left_join_table_1) );
                vec_1__df3_y->push_back( _df1_y.ARRAYELEM(tmp_left_join_table_1) );
                vec_1__df3_x1->push_back( _df2_x1.ARRAYELEM(right_join_table_1) );
                vec_1__df3_y1->push_back( _df2_y1.ARRAYELEM(right_join_table_1) );
                tmp_left_join_table_1++;
                table_new_counter_join_1++;
            }
            int tmp_right_join_table_1 = right_join_table_1 + 1 ;
            while((tmp_right_join_table_1 < rsize_t2_1 + 1) && (_df1_id.ARRAYELEM(left_join_table_1) == _df2_id.ARRAYELEM(tmp_right_join_table_1)))
            {
                vec_1__df3_id->push_back( _df1_id.ARRAYELEM(left_join_table_1) );
                vec_1__df3_x->push_back( _df1_x.ARRAYELEM(left_join_table_1) );
                vec_1__df3_y->push_back( _df1_y.ARRAYELEM(left_join_table_1) );
                vec_1__df3_x1->push_back( _df2_x1.ARRAYELEM(tmp_right_join_table_1) );
                vec_1__df3_y1->push_back( _df2_y1.ARRAYELEM(tmp_right_join_table_1) );
                tmp_right_join_table_1++;
                table_new_counter_join_1++;
            }
            left_join_table_1++;
            right_join_table_1++;
        }
        else if (_df1_id.ARRAYELEM(left_join_table_1) < _df2_id.ARRAYELEM(right_join_table_1))
            left_join_table_1++;
        else
            right_join_table_1++;
    }
    _df3_id = j2c_array<int64_t>::new_j2c_array_1d(vec_1__df3_id->data(), vec_1__df3_id->size() );
    _df3_x = j2c_array<double>::new_j2c_array_1d(vec_1__df3_x->data(), vec_1__df3_x->size() );
    _df3_y = j2c_array<double>::new_j2c_array_1d(vec_1__df3_y->data(), vec_1__df3_y->size() );
    _df3_x1 = j2c_array<double>::new_j2c_array_1d(vec_1__df3_x1->data(), vec_1__df3_x1->size() );
    _df3_y1 = j2c_array<double>::new_j2c_array_1d(vec_1__df3_y1->data(), vec_1__df3_y1->size() );
    ;
    t2 = MPI_Wtime();
    parallel_ir_save_array_len_1_1 = _df3_x.ARRAYSIZE(1);
    parallel_ir_reduction_output_1 = 0.0;
    for ( parfor_index_1_1 = 1; parfor_index_1_1 <= (int64_t)parallel_ir_save_array_len_1_1; parfor_index_1_1 += 1)
    {
        ;
        parallel_ir_array_temp__20_3_1 = _df3_x.ARRAYELEM(parfor_index_1_1);
        SSAValue3 = (parallel_ir_reduction_output_1) + (parallel_ir_array_temp__20_3_1);
        parallel_ir_reduction_output_1 = SSAValue3;
    }
    ;
    SSAValue2 = parallel_ir_reduction_output_1;
    *ret0 = SSAValue2;
    return;

}


void ppcumsum_testp271_unaliased(ASCIIString& __restrict file_name, ASCIIString& __restrict file_name2, double * __restrict ret0)
{
    pppcumsum_testp271 pselfp;
    j2c_array< int64_t >  _df1_id;
    j2c_array< double >  _df1_x;
    j2c_array< double >  _df1_y;
    j2c_array< int64_t >  _df2_id;
    j2c_array< double >  _df2_x1;
    j2c_array< double >  _df2_y1;
    j2c_array< int64_t >  _df3_id;
    j2c_array< double >  _df3_x;
    j2c_array< double >  _df3_y;
    j2c_array< double >  _df3_x1;
    j2c_array< double >  _df3_y1;
    int64_t __hpat_h5_dim_size_2_1;
    int64_t __hpat_h5_dim_size_5_1;
    hsize_t* SSAValue0;
    hsize_t* SSAValue1;
    double SSAValue2;
    int64_t parfor_index_1_1;
    double parallel_ir_array_temp__20_3_1;
    int64_t parallel_ir_save_array_len_1_1;
    double parallel_ir_reduction_output_1;
    double SSAValue3;
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
    int64_t __hpat_dist_arr_start_5;
    int64_t __hpat_dist_arr_div_5;
    int64_t __hpat_dist_arr_count_5;
    int64_t __hpat_dist_arr_start_6;
    int64_t __hpat_dist_arr_div_6;
    int64_t __hpat_dist_arr_count_6;
    ;;
    MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
    MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
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
    dataset_id_2 = H5Dopen2(file_id_2, "/id", H5P_DEFAULT);
    assert(dataset_id_2 != -1);
    ;
    hid_t space_id_2 = H5Dget_space(dataset_id_2);
    assert(space_id_2 != -1);
    hsize_t data_ndim_2 = H5Sget_simple_extent_ndims(space_id_2);
    hsize_t space_dims_2[data_ndim_2];
    H5Sget_simple_extent_dims(space_id_2, space_dims_2, NULL);
    SSAValue0 = space_dims_2;;
    __hpat_h5_dim_size_2_1 = SSAValue0[1-1];
    __hpat_dist_arr_div_1 = (__hpat_h5_dim_size_2_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_1 = (__hpat_node_id) * (__hpat_dist_arr_div_1);
    __hpat_dist_arr_count_1 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_2_1-__hpat_node_id*__hpat_dist_arr_div_1 : __hpat_dist_arr_div_1);
    _df1_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
    hsize_t CGen_HDF5_start_2[data_ndim_2];
    hsize_t CGen_HDF5_count_2[data_ndim_2];
    CGen_HDF5_start_2[0] = __hpat_dist_arr_start_1;
    CGen_HDF5_count_2[0] = __hpat_dist_arr_count_1;
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
    ret_2 = H5Dread(dataset_id_2, H5T_NATIVE_LLONG, mem_dataspace_2, space_id_2, xfer_plist_2, _df1_id.getData());
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
    dataset_id_3 = H5Dopen2(file_id_3, "/x", H5P_DEFAULT);
    assert(dataset_id_3 != -1);
    ;
    hid_t space_id_3 = H5Dget_space(dataset_id_3);
    assert(space_id_3 != -1);
    hsize_t data_ndim_3 = H5Sget_simple_extent_ndims(space_id_3);
    hsize_t space_dims_3[data_ndim_3];
    H5Sget_simple_extent_dims(space_id_3, space_dims_3, NULL);
    ;
    __hpat_dist_arr_div_2 = (__hpat_h5_dim_size_2_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_2 = (__hpat_node_id) * (__hpat_dist_arr_div_2);
    __hpat_dist_arr_count_2 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_2_1-__hpat_node_id*__hpat_dist_arr_div_2 : __hpat_dist_arr_div_2);
    _df1_x = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_2);
    hsize_t CGen_HDF5_start_3[data_ndim_3];
    hsize_t CGen_HDF5_count_3[data_ndim_3];
    CGen_HDF5_start_3[0] = __hpat_dist_arr_start_2;
    CGen_HDF5_count_3[0] = __hpat_dist_arr_count_2;
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
    ret_3 = H5Dread(dataset_id_3, H5T_NATIVE_DOUBLE, mem_dataspace_3, space_id_3, xfer_plist_3, _df1_x.getData());
    assert(ret_3 != -1);
    ;
    ;
    H5Dclose(dataset_id_3);
    H5Fclose(file_id_3);
    ;
    hid_t plist_id_4 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_4 != -1);
    herr_t ret_4;
    hid_t file_id_4;
    ret_4 = H5Pset_fapl_mpio(plist_id_4, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_4 != -1);
    file_id_4 = H5Fopen((const char*)file_name.data.data, H5F_ACC_RDONLY, plist_id_4);
    assert(file_id_4 != -1);
    ret_4 = H5Pclose(plist_id_4);
    assert(ret_4 != -1);
    hid_t dataset_id_4;
    dataset_id_4 = H5Dopen2(file_id_4, "/y", H5P_DEFAULT);
    assert(dataset_id_4 != -1);
    ;
    hid_t space_id_4 = H5Dget_space(dataset_id_4);
    assert(space_id_4 != -1);
    hsize_t data_ndim_4 = H5Sget_simple_extent_ndims(space_id_4);
    hsize_t space_dims_4[data_ndim_4];
    H5Sget_simple_extent_dims(space_id_4, space_dims_4, NULL);
    ;
    __hpat_dist_arr_div_3 = (__hpat_h5_dim_size_2_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_3 = (__hpat_node_id) * (__hpat_dist_arr_div_3);
    __hpat_dist_arr_count_3 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_2_1-__hpat_node_id*__hpat_dist_arr_div_3 : __hpat_dist_arr_div_3);
    _df1_y = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_3);
    hsize_t CGen_HDF5_start_4[data_ndim_4];
    hsize_t CGen_HDF5_count_4[data_ndim_4];
    CGen_HDF5_start_4[0] = __hpat_dist_arr_start_3;
    CGen_HDF5_count_4[0] = __hpat_dist_arr_count_3;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_4; i_CGen_dim++)
    {
        CGen_HDF5_start_4[i_CGen_dim] = 0;
        CGen_HDF5_count_4[i_CGen_dim] = space_dims_4[i_CGen_dim];
    }
    ret_4 = H5Sselect_hyperslab(space_id_4, H5S_SELECT_SET, CGen_HDF5_start_4, NULL, CGen_HDF5_count_4, NULL);
    assert(ret_4 != -1);
    hid_t mem_dataspace_4 = H5Screate_simple (data_ndim_4, CGen_HDF5_count_4, NULL);
    assert (mem_dataspace_4 != -1);
    hid_t xfer_plist_4 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_4 != -1);
    double h5_read_start_4 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_4, H5FD_MPIO_COLLECTIVE);
    ret_4 = H5Dread(dataset_id_4, H5T_NATIVE_DOUBLE, mem_dataspace_4, space_id_4, xfer_plist_4, _df1_y.getData());
    assert(ret_4 != -1);
    ;
    ;
    H5Dclose(dataset_id_4);
    H5Fclose(file_id_4);
    ;
    hid_t plist_id_5 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_5 != -1);
    herr_t ret_5;
    hid_t file_id_5;
    ret_5 = H5Pset_fapl_mpio(plist_id_5, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_5 != -1);
    file_id_5 = H5Fopen((const char*)file_name2.data.data, H5F_ACC_RDONLY, plist_id_5);
    assert(file_id_5 != -1);
    ret_5 = H5Pclose(plist_id_5);
    assert(ret_5 != -1);
    hid_t dataset_id_5;
    dataset_id_5 = H5Dopen2(file_id_5, "/id", H5P_DEFAULT);
    assert(dataset_id_5 != -1);
    ;
    hid_t space_id_5 = H5Dget_space(dataset_id_5);
    assert(space_id_5 != -1);
    hsize_t data_ndim_5 = H5Sget_simple_extent_ndims(space_id_5);
    hsize_t space_dims_5[data_ndim_5];
    H5Sget_simple_extent_dims(space_id_5, space_dims_5, NULL);
    SSAValue1 = space_dims_5;;
    __hpat_h5_dim_size_5_1 = SSAValue1[1-1];
    __hpat_dist_arr_div_4 = (__hpat_h5_dim_size_5_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_4 = (__hpat_node_id) * (__hpat_dist_arr_div_4);
    __hpat_dist_arr_count_4 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_5_1-__hpat_node_id*__hpat_dist_arr_div_4 : __hpat_dist_arr_div_4);
    _df2_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_4);
    hsize_t CGen_HDF5_start_5[data_ndim_5];
    hsize_t CGen_HDF5_count_5[data_ndim_5];
    CGen_HDF5_start_5[0] = __hpat_dist_arr_start_4;
    CGen_HDF5_count_5[0] = __hpat_dist_arr_count_4;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_5; i_CGen_dim++)
    {
        CGen_HDF5_start_5[i_CGen_dim] = 0;
        CGen_HDF5_count_5[i_CGen_dim] = space_dims_5[i_CGen_dim];
    }
    ret_5 = H5Sselect_hyperslab(space_id_5, H5S_SELECT_SET, CGen_HDF5_start_5, NULL, CGen_HDF5_count_5, NULL);
    assert(ret_5 != -1);
    hid_t mem_dataspace_5 = H5Screate_simple (data_ndim_5, CGen_HDF5_count_5, NULL);
    assert (mem_dataspace_5 != -1);
    hid_t xfer_plist_5 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_5 != -1);
    double h5_read_start_5 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_5, H5FD_MPIO_COLLECTIVE);
    ret_5 = H5Dread(dataset_id_5, H5T_NATIVE_LLONG, mem_dataspace_5, space_id_5, xfer_plist_5, _df2_id.getData());
    assert(ret_5 != -1);
    ;
    ;
    H5Dclose(dataset_id_5);
    H5Fclose(file_id_5);
    ;
    hid_t plist_id_6 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_6 != -1);
    herr_t ret_6;
    hid_t file_id_6;
    ret_6 = H5Pset_fapl_mpio(plist_id_6, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_6 != -1);
    file_id_6 = H5Fopen((const char*)file_name2.data.data, H5F_ACC_RDONLY, plist_id_6);
    assert(file_id_6 != -1);
    ret_6 = H5Pclose(plist_id_6);
    assert(ret_6 != -1);
    hid_t dataset_id_6;
    dataset_id_6 = H5Dopen2(file_id_6, "/x1", H5P_DEFAULT);
    assert(dataset_id_6 != -1);
    ;
    hid_t space_id_6 = H5Dget_space(dataset_id_6);
    assert(space_id_6 != -1);
    hsize_t data_ndim_6 = H5Sget_simple_extent_ndims(space_id_6);
    hsize_t space_dims_6[data_ndim_6];
    H5Sget_simple_extent_dims(space_id_6, space_dims_6, NULL);
    ;
    __hpat_dist_arr_div_5 = (__hpat_h5_dim_size_5_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_5 = (__hpat_node_id) * (__hpat_dist_arr_div_5);
    __hpat_dist_arr_count_5 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_5_1-__hpat_node_id*__hpat_dist_arr_div_5 : __hpat_dist_arr_div_5);
    _df2_x1 = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_5);
    hsize_t CGen_HDF5_start_6[data_ndim_6];
    hsize_t CGen_HDF5_count_6[data_ndim_6];
    CGen_HDF5_start_6[0] = __hpat_dist_arr_start_5;
    CGen_HDF5_count_6[0] = __hpat_dist_arr_count_5;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_6; i_CGen_dim++)
    {
        CGen_HDF5_start_6[i_CGen_dim] = 0;
        CGen_HDF5_count_6[i_CGen_dim] = space_dims_6[i_CGen_dim];
    }
    ret_6 = H5Sselect_hyperslab(space_id_6, H5S_SELECT_SET, CGen_HDF5_start_6, NULL, CGen_HDF5_count_6, NULL);
    assert(ret_6 != -1);
    hid_t mem_dataspace_6 = H5Screate_simple (data_ndim_6, CGen_HDF5_count_6, NULL);
    assert (mem_dataspace_6 != -1);
    hid_t xfer_plist_6 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_6 != -1);
    double h5_read_start_6 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_6, H5FD_MPIO_COLLECTIVE);
    ret_6 = H5Dread(dataset_id_6, H5T_NATIVE_DOUBLE, mem_dataspace_6, space_id_6, xfer_plist_6, _df2_x1.getData());
    assert(ret_6 != -1);
    ;
    ;
    H5Dclose(dataset_id_6);
    H5Fclose(file_id_6);
    ;
    hid_t plist_id_7 = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist_id_7 != -1);
    herr_t ret_7;
    hid_t file_id_7;
    ret_7 = H5Pset_fapl_mpio(plist_id_7, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(ret_7 != -1);
    file_id_7 = H5Fopen((const char*)file_name2.data.data, H5F_ACC_RDONLY, plist_id_7);
    assert(file_id_7 != -1);
    ret_7 = H5Pclose(plist_id_7);
    assert(ret_7 != -1);
    hid_t dataset_id_7;
    dataset_id_7 = H5Dopen2(file_id_7, "/y1", H5P_DEFAULT);
    assert(dataset_id_7 != -1);
    ;
    hid_t space_id_7 = H5Dget_space(dataset_id_7);
    assert(space_id_7 != -1);
    hsize_t data_ndim_7 = H5Sget_simple_extent_ndims(space_id_7);
    hsize_t space_dims_7[data_ndim_7];
    H5Sget_simple_extent_dims(space_id_7, space_dims_7, NULL);
    ;
    __hpat_dist_arr_div_6 = (__hpat_h5_dim_size_5_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_6 = (__hpat_node_id) * (__hpat_dist_arr_div_6);
    __hpat_dist_arr_count_6 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_5_1-__hpat_node_id*__hpat_dist_arr_div_6 : __hpat_dist_arr_div_6);
    _df2_y1 = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_6);
    hsize_t CGen_HDF5_start_7[data_ndim_7];
    hsize_t CGen_HDF5_count_7[data_ndim_7];
    CGen_HDF5_start_7[0] = __hpat_dist_arr_start_6;
    CGen_HDF5_count_7[0] = __hpat_dist_arr_count_6;
    for(int i_CGen_dim=1; i_CGen_dim<data_ndim_7; i_CGen_dim++)
    {
        CGen_HDF5_start_7[i_CGen_dim] = 0;
        CGen_HDF5_count_7[i_CGen_dim] = space_dims_7[i_CGen_dim];
    }
    ret_7 = H5Sselect_hyperslab(space_id_7, H5S_SELECT_SET, CGen_HDF5_start_7, NULL, CGen_HDF5_count_7, NULL);
    assert(ret_7 != -1);
    hid_t mem_dataspace_7 = H5Screate_simple (data_ndim_7, CGen_HDF5_count_7, NULL);
    assert (mem_dataspace_7 != -1);
    hid_t xfer_plist_7 = H5Pcreate (H5P_DATASET_XFER);
    assert(xfer_plist_7 != -1);
    double h5_read_start_7 = MPI_Wtime();
    H5Pset_dxpl_mpio(xfer_plist_7, H5FD_MPIO_COLLECTIVE);
    ret_7 = H5Dread(dataset_id_7, H5T_NATIVE_DOUBLE, mem_dataspace_7, space_id_7, xfer_plist_7, _df2_y1.getData());
    assert(ret_7 != -1);
    ;
    ;
    H5Dclose(dataset_id_7);
    H5Fclose(file_id_7);
    ;
    int join_num_pes_1;
    MPI_Comm join_comm_1;
    join_num_pes_1 = __hpat_num_pes ;
    join_comm_1 = MPI_COMM_WORLD ;
    int * scount_t1_1;
    int * scount_t2_1;
    int * scount_t1_tmp_1;
    int * scount_t2_tmp_1;
    scount_t1_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    memset (scount_t1_1, 0, sizeof(int)* join_num_pes_1);
    scount_t1_tmp_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    memset (scount_t1_tmp_1, 0, sizeof(int)* join_num_pes_1);
    scount_t2_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    memset (scount_t2_1, 0, sizeof(int)* join_num_pes_1);
    scount_t2_tmp_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    memset (scount_t2_tmp_1, 0, sizeof(int)* join_num_pes_1);
    int rsize_t1_1 = 0;
    int rsize_t2_1 = 0;
    int * rcount_t1_1;
    int * rcount_t2_1;
    rcount_t1_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    rcount_t2_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    int * sdis_t1_1;
    int * rdis_t1_1;
    sdis_t1_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    rdis_t1_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    int *sdis_t2_1;
    int *rdis_t2_1;
    sdis_t2_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    rdis_t2_1 = (int*)malloc(sizeof(int)* join_num_pes_1);
    int t1c1_length_join_1 = _df1_id.ARRAYLEN() ;
    int t2c1_length_join_1 = _df2_id.ARRAYLEN() ;
    for (int i = 1 ; i <  t1c1_length_join_1 + 1 ; i++)
    {
        int node_id = _df1_id.ARRAYELEM(i) % join_num_pes_1 ;
        scount_t1_1[node_id]++;
    }
    sdis_t1_1[0]=0;
    for(int i=1;i < join_num_pes_1 ;i++)
    {
        sdis_t1_1[i]=scount_t1_1[i-1] + sdis_t1_1[i-1];
    }
    MPI_Alltoall(scount_t1_1,1,MPI_INT,rcount_t1_1,1,MPI_INT, join_comm_1);
    j2c_array< int64_t > _df1_id_tmp_join_1 = j2c_array<int64_t>::new_j2c_array_1d(NULL, t1c1_length_join_1 );
    j2c_array< double > _df1_x_tmp_join_1 = j2c_array<double>::new_j2c_array_1d(NULL, t1c1_length_join_1 );
    j2c_array< double > _df1_y_tmp_join_1 = j2c_array<double>::new_j2c_array_1d(NULL, t1c1_length_join_1 );
    for (int i = 1 ; i <  t1c1_length_join_1 + 1 ; i++)
    {
        int node_id = _df1_id.ARRAYELEM(i) % join_num_pes_1;
        _df1_id_tmp_join_1.ARRAYELEM(sdis_t1_1[node_id]+scount_t1_tmp_1[node_id]+1) = _df1_id.ARRAYELEM(i);
        _df1_x_tmp_join_1.ARRAYELEM(sdis_t1_1[node_id]+scount_t1_tmp_1[node_id]+1) = _df1_x.ARRAYELEM(i);
        _df1_y_tmp_join_1.ARRAYELEM(sdis_t1_1[node_id]+scount_t1_tmp_1[node_id]+1) = _df1_y.ARRAYELEM(i);
        scount_t1_tmp_1[node_id]++;
    }
    for (int i = 1 ; i <  t2c1_length_join_1 + 1 ; i++)
    {
        int node_id = _df2_id.ARRAYELEM(i) % join_num_pes_1 ;
        scount_t2_1[node_id]++;
    }
    sdis_t2_1[0]=0;
    for(int i=1;i < __hpat_num_pes;i++)
    {
        sdis_t2_1[i]=scount_t2_1[i-1] + sdis_t2_1[i-1];
    }
    MPI_Alltoall(scount_t2_1,1,MPI_INT, rcount_t2_1,1,MPI_INT,MPI_COMM_WORLD);
    j2c_array< int64_t > _df2_id_tmp_join_1 = j2c_array<int64_t>::new_j2c_array_1d(NULL, t2c1_length_join_1);
    j2c_array< double > _df2_x1_tmp_join_1 = j2c_array<double>::new_j2c_array_1d(NULL, t2c1_length_join_1);
    j2c_array< double > _df2_y1_tmp_join_1 = j2c_array<double>::new_j2c_array_1d(NULL, t2c1_length_join_1);
    for (int i = 1 ; i <   t2c1_length_join_1 + 1 ; i++)
    {
        int node_id = _df2_id.ARRAYELEM(i) % join_num_pes_1 ;
        _df2_id_tmp_join_1.ARRAYELEM(sdis_t2_1[node_id]+scount_t2_tmp_1[node_id]+1) = _df2_id.ARRAYELEM(i);
        _df2_x1_tmp_join_1.ARRAYELEM(sdis_t2_1[node_id]+scount_t2_tmp_1[node_id]+1) = _df2_x1.ARRAYELEM(i);
        _df2_y1_tmp_join_1.ARRAYELEM(sdis_t2_1[node_id]+scount_t2_tmp_1[node_id]+1) = _df2_y1.ARRAYELEM(i);
        scount_t2_tmp_1[node_id]++;
    }
    rdis_t1_1[0]=0;
    rdis_t2_1[0]=0;
    for(int i=1;i < join_num_pes_1 ;i++)
    {
        rdis_t1_1[i] = rcount_t1_1[i-1] + rdis_t1_1[i-1];
        rdis_t2_1[i] = rcount_t2_1[i-1] + rdis_t2_1[i-1];
    }
    for(int i=0;i< join_num_pes_1 ;i++)
    {
        rsize_t1_1 = rsize_t1_1 + rcount_t1_1[i];
        rsize_t2_1 = rsize_t2_1 + rcount_t2_1[i];
    }
    j2c_array< int64_t > rbuf__df1_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, rsize_t1_1);
    MPI_Alltoallv(_df1_id_tmp_join_1.getData(), scount_t1_1, sdis_t1_1, MPI_INT64_T,
        rbuf__df1_id.getData(), rcount_t1_1, rdis_t1_1, MPI_INT64_T, join_comm_1);
    _df1_id = rbuf__df1_id;
    j2c_array< double > rbuf__df1_x = j2c_array<double>::new_j2c_array_1d(NULL, rsize_t1_1);
    MPI_Alltoallv(_df1_x_tmp_join_1.getData(), scount_t1_1, sdis_t1_1, MPI_INT64_T,
        rbuf__df1_x.getData(), rcount_t1_1, rdis_t1_1, MPI_INT64_T, join_comm_1);
    _df1_x = rbuf__df1_x;
    j2c_array< double > rbuf__df1_y = j2c_array<double>::new_j2c_array_1d(NULL, rsize_t1_1);
    MPI_Alltoallv(_df1_y_tmp_join_1.getData(), scount_t1_1, sdis_t1_1, MPI_INT64_T,
        rbuf__df1_y.getData(), rcount_t1_1, rdis_t1_1, MPI_INT64_T, join_comm_1);
    _df1_y = rbuf__df1_y;
    j2c_array< int64_t > rbuf__df2_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, rsize_t2_1);
    MPI_Alltoallv(_df2_id_tmp_join_1.getData(), scount_t2_1, sdis_t2_1, MPI_INT64_T,
        rbuf__df2_id.getData(), rcount_t2_1, rdis_t2_1, MPI_INT64_T, join_comm_1);
    _df2_id = rbuf__df2_id;
    j2c_array< double > rbuf__df2_x1 = j2c_array<double>::new_j2c_array_1d(NULL, rsize_t2_1);
    MPI_Alltoallv(_df2_x1_tmp_join_1.getData(), scount_t2_1, sdis_t2_1, MPI_INT64_T,
        rbuf__df2_x1.getData(), rcount_t2_1, rdis_t2_1, MPI_INT64_T, join_comm_1);
    _df2_x1 = rbuf__df2_x1;
    j2c_array< double > rbuf__df2_y1 = j2c_array<double>::new_j2c_array_1d(NULL, rsize_t2_1);
    MPI_Alltoallv(_df2_y1_tmp_join_1.getData(), scount_t2_1, sdis_t2_1, MPI_INT64_T,
        rbuf__df2_y1.getData(), rcount_t2_1, rdis_t2_1, MPI_INT64_T, join_comm_1);
    _df2_y1 = rbuf__df2_y1;
    int table_new_counter_join_1 = 1 ;
    std::vector< int64_t> *vec_1__df3_id = new std::vector< int64_t>();
    std::vector< double> *vec_1__df3_x = new std::vector< double>();
    std::vector< double> *vec_1__df3_y = new std::vector< double>();
    std::vector< double> *vec_1__df3_x1 = new std::vector< double>();
    std::vector< double> *vec_1__df3_y1 = new std::vector< double>();
    int64_t * t1_all_arrays_1[3 - 1];
    int64_t * t2_all_arrays_1[3 - 1];
    t1_all_arrays_1 [0] = ( int64_t *) _df1_x.getData();
    t1_all_arrays_1 [1] = ( int64_t *) _df1_y.getData();
    t2_all_arrays_1 [0] = ( int64_t *) _df2_x1.getData();
    t2_all_arrays_1 [1] = ( int64_t *) _df2_y1.getData();
    __hpat_timsort(( int64_t *) _df1_id.getData(), rsize_t1_1 , t1_all_arrays_1, 3 - 1);
    __hpat_timsort(( int64_t *) _df2_id.getData(), rsize_t2_1 , t2_all_arrays_1, 3 - 1);
    int left_join_table_1 = 1;
    int right_join_table_1 = 1;
    while ( (left_join_table_1 < rsize_t1_1 + 1) && (right_join_table_1 < rsize_t2_1 + 1) )
    {
        if(_df1_id.ARRAYELEM(left_join_table_1) == _df2_id.ARRAYELEM(right_join_table_1))
        {
            vec_1__df3_id->push_back( _df1_id.ARRAYELEM(left_join_table_1) );
            vec_1__df3_x->push_back( _df1_x.ARRAYELEM(left_join_table_1) );
            vec_1__df3_y->push_back( _df1_y.ARRAYELEM(left_join_table_1) );
            vec_1__df3_x1->push_back( _df2_x1.ARRAYELEM(right_join_table_1) );
            vec_1__df3_y1->push_back( _df2_y1.ARRAYELEM(right_join_table_1) );
            table_new_counter_join_1++;
            int tmp_left_join_table_1 = left_join_table_1 + 1 ;
            while((tmp_left_join_table_1 < rsize_t1_1 + 1) && (_df1_id.ARRAYELEM(tmp_left_join_table_1) == _df2_id.ARRAYELEM(right_join_table_1)))
            {
                vec_1__df3_id->push_back( _df1_id.ARRAYELEM(tmp_left_join_table_1) );
                vec_1__df3_x->push_back( _df1_x.ARRAYELEM(tmp_left_join_table_1) );
                vec_1__df3_y->push_back( _df1_y.ARRAYELEM(tmp_left_join_table_1) );
                vec_1__df3_x1->push_back( _df2_x1.ARRAYELEM(right_join_table_1) );
                vec_1__df3_y1->push_back( _df2_y1.ARRAYELEM(right_join_table_1) );
                tmp_left_join_table_1++;
                table_new_counter_join_1++;
            }
            int tmp_right_join_table_1 = right_join_table_1 + 1 ;
            while((tmp_right_join_table_1 < rsize_t2_1 + 1) && (_df1_id.ARRAYELEM(left_join_table_1) == _df2_id.ARRAYELEM(tmp_right_join_table_1)))
            {
                vec_1__df3_id->push_back( _df1_id.ARRAYELEM(left_join_table_1) );
                vec_1__df3_x->push_back( _df1_x.ARRAYELEM(left_join_table_1) );
                vec_1__df3_y->push_back( _df1_y.ARRAYELEM(left_join_table_1) );
                vec_1__df3_x1->push_back( _df2_x1.ARRAYELEM(tmp_right_join_table_1) );
                vec_1__df3_y1->push_back( _df2_y1.ARRAYELEM(tmp_right_join_table_1) );
                tmp_right_join_table_1++;
                table_new_counter_join_1++;
            }
            left_join_table_1++;
            right_join_table_1++;
        }
        else if (_df1_id.ARRAYELEM(left_join_table_1) < _df2_id.ARRAYELEM(right_join_table_1))
            left_join_table_1++;
        else
            right_join_table_1++;
    }
    _df3_id = j2c_array<int64_t>::new_j2c_array_1d(vec_1__df3_id->data(), vec_1__df3_id->size() );
    _df3_x = j2c_array<double>::new_j2c_array_1d(vec_1__df3_x->data(), vec_1__df3_x->size() );
    _df3_y = j2c_array<double>::new_j2c_array_1d(vec_1__df3_y->data(), vec_1__df3_y->size() );
    _df3_x1 = j2c_array<double>::new_j2c_array_1d(vec_1__df3_x1->data(), vec_1__df3_x1->size() );
    _df3_y1 = j2c_array<double>::new_j2c_array_1d(vec_1__df3_y1->data(), vec_1__df3_y1->size() );
    ;
    parallel_ir_save_array_len_1_1 = _df3_x.ARRAYSIZE(1);
    parallel_ir_reduction_output_1 = 0.0;
    for ( parfor_index_1_1 = 1; parfor_index_1_1 <= (int64_t)parallel_ir_save_array_len_1_1; parfor_index_1_1 += 1)
    {
        ;
        parallel_ir_array_temp__20_3_1 = _df3_x.ARRAYELEM(parfor_index_1_1);
        SSAValue3 = (parallel_ir_reduction_output_1) + (parallel_ir_array_temp__20_3_1);
        parallel_ir_reduction_output_1 = SSAValue3;
    }
    ;
    SSAValue2 = parallel_ir_reduction_output_1;
    *ret0 = SSAValue2;
    return;

}


extern "C" void _ppcumsum_testp271_unaliased_(int run_where, ASCIIString& __restrict file_name, ASCIIString& __restrict file_name2 , double* __restrict ret0 , bool genMain = true)
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
        mainFileData << file_name2 << std::endl;
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
        mainFile << "    ASCIIString file_name2;" << std::endl;
        mainFile << "    mainFileData >> file_name2;" << std::endl;
        mainFile << "    _ppcumsum_testp271_unaliased_(runwhere, file_name, file_name2, &ret0, false);" << std::endl;
        mainFile << "    MPI_Finalize();" << std::endl;
        mainFile << "    return 0;" << std::endl;
        mainFile << "}" << std::endl;
        mainFile.close();
        std::ofstream mainFileSh(newMainSh.str());
        mainFileSh << "#!/bin/sh" << std::endl;
        mainFileSh << "mpiicpc -O3 -std=c++11  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << " -lhdf5  -lm " << std::endl;
        mainFileSh.close();
    }

    ppcumsum_testp271_unaliased(file_name, file_name2, ret0);
}


extern "C" void _ppcumsum_testp271_(int run_where, ASCIIString&  file_name, ASCIIString&  file_name2 , double* __restrict ret0 , bool genMain = true)
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
        mainFileData << file_name2 << std::endl;
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
        mainFile << "    ASCIIString file_name2;" << std::endl;
        mainFile << "    mainFileData >> file_name2;" << std::endl;
        mainFile << "    _ppcumsum_testp271_(runwhere, file_name, file_name2, &ret0, false);" << std::endl;
        mainFile << "    MPI_Finalize();" << std::endl;
        mainFile << "    return 0;" << std::endl;
        mainFile << "}" << std::endl;
        mainFile.close();
        std::ofstream mainFileSh(newMainSh.str());
        mainFileSh << "#!/bin/sh" << std::endl;
        mainFileSh << "mpiicpc -O3 -std=c++11  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << " -lhdf5  -lm " << std::endl;
        mainFileSh.close();
    }

    ppcumsum_testp271(file_name, file_name2, ret0);
}


extern "C"
void *j2c_array_new(int key, void*data, unsigned ndim, int64_t *dims)
{
    void *a = NULL;
    switch(key)
    {
        case 2:
            a = new  j2c_array< uint8_t > ((uint8_t*)data, ndim, dims);
            break;
        default:
            fprintf(stderr, "j2c_array_new called with invalid key %d", key);
            assert(false);
            break;
    }
    return a;
}
