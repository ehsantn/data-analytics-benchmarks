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
    j2c_array< bool >  _df_cond_e;
    j2c_array< int64_t >  _df1_id;
    j2c_array< double >  _df1_x;
    j2c_array< double >  _df1_y;
    int64_t __hpat_h5_dim_size_2_1;
    hsize_t* SSAValue0;
    double SSAValue1;
    int64_t parallel_ir_array_temp__3_4_1;
    int64_t parfor_index_1_3;
    int64_t parallel_ir_save_array_len_1_3;
    bool SSAValue2;
    j2c_array< bool >  parallel_ir_new_array_name_3_1;
    bool parallel_ir_array_temp__24_6_1;
    int64_t parfor_index_1_7;
    double parallel_ir_array_temp__13_9_1;
    int64_t parallel_ir_save_array_len_1_7;
    double parallel_ir_reduction_output_7;
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
    int64_t __hpat_loop_start_3;
    int64_t __hpat_loop_end_3;
    int64_t __hpat_loop_div_3;
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
    _df_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
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
    ret_2 = H5Dread(dataset_id_2, H5T_NATIVE_LLONG, mem_dataspace_2, space_id_2, xfer_plist_2, _df_id.getData());
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
    _df_x = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_2);
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
    ret_3 = H5Dread(dataset_id_3, H5T_NATIVE_DOUBLE, mem_dataspace_3, space_id_3, xfer_plist_3, _df_x.getData());
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
    _df_y = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_3);
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
    ret_4 = H5Dread(dataset_id_4, H5T_NATIVE_DOUBLE, mem_dataspace_4, space_id_4, xfer_plist_4, _df_y.getData());
    assert(ret_4 != -1);
    ;
    ;
    H5Dclose(dataset_id_4);
    H5Fclose(file_id_4);
    ;
    t1 = MPI_Wtime();
    parallel_ir_save_array_len_1_3 = __hpat_h5_dim_size_2_1;
    __hpat_dist_arr_div_4 = (parallel_ir_save_array_len_1_3) / (__hpat_num_pes);
    __hpat_dist_arr_start_4 = (__hpat_node_id) * (__hpat_dist_arr_div_4);
    __hpat_dist_arr_count_4 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_3-__hpat_node_id*__hpat_dist_arr_div_4 : __hpat_dist_arr_div_4);
    parallel_ir_new_array_name_3_1 = j2c_array<bool>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_4);
    __hpat_loop_div_3 = (parallel_ir_save_array_len_1_3) / (__hpat_num_pes);
    __hpat_loop_start_3 = ((__hpat_node_id) * (__hpat_loop_div_3)) + (1);
    __hpat_loop_end_3 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_3 : (__hpat_node_id+1)*__hpat_loop_div_3);
    for ( parfor_index_1_3 = __hpat_loop_start_3; parfor_index_1_3 <= (int64_t)__hpat_loop_end_3; parfor_index_1_3 += 1)
    {
        ;
        parallel_ir_array_temp__3_4_1 = _df_id.ARRAYELEM(((parfor_index_1_3) - (__hpat_loop_start_3)) + (1));
        SSAValue2 = (parallel_ir_array_temp__3_4_1) < (100);
        parallel_ir_array_temp__24_6_1 = SSAValue2;
        parallel_ir_new_array_name_3_1.ARRAYELEM(((parfor_index_1_3) - (__hpat_loop_start_3)) + (1)) = parallel_ir_array_temp__24_6_1;
    }
    ;
    _df_cond_e = parallel_ir_new_array_name_3_1;
    int _df_id_array_length_filter1 = _df_id.ARRAYLEN();
    int write_index_filter1 = 1;
    for (int index = 1 ; index < _df_id_array_length_filter1 + 1 ; index++)
    {
        if ( _df_cond_e.ARRAYELEM(index) )
        {
            _df_id.ARRAYELEM(write_index_filter1) =  _df_id.ARRAYELEM(index);
            _df_x.ARRAYELEM(write_index_filter1) =  _df_x.ARRAYELEM(index);
            _df_y.ARRAYELEM(write_index_filter1) =  _df_y.ARRAYELEM(index);
            write_index_filter1 = write_index_filter1 + 1;
        };
    };
    _df_id.dims[0] =  write_index_filter1 - 1;
    _df1_id = _df_id ;
    _df_x.dims[0] =  write_index_filter1 - 1;
    _df1_x = _df_x ;
    _df_y.dims[0] =  write_index_filter1 - 1;
    _df1_y = _df_y ;
    t2 = MPI_Wtime();
    ;
    parallel_ir_save_array_len_1_7 = _df1_x.ARRAYSIZE(1);
    parallel_ir_reduction_output_7 = 0.0;
    for ( parfor_index_1_7 = 1; parfor_index_1_7 <= (int64_t)parallel_ir_save_array_len_1_7; parfor_index_1_7 += 1)
    {
        ;
        parallel_ir_array_temp__13_9_1 = _df1_x.ARRAYELEM(parfor_index_1_7);
        SSAValue3 = (parallel_ir_reduction_output_7) + (parallel_ir_array_temp__13_9_1);
        parallel_ir_reduction_output_7 = SSAValue3;
    }
    ;
    SSAValue1 = parallel_ir_reduction_output_7;
    *ret0 = SSAValue1;
    return;

}


void ppcumsum_testp271_unaliased(ASCIIString& __restrict file_name, double * __restrict ret0)
{
    pppcumsum_testp271 pselfp;
    j2c_array< int64_t >  _df_id;
    j2c_array< double >  _df_x;
    j2c_array< double >  _df_y;
    j2c_array< bool >  _df_cond_e;
    j2c_array< int64_t >  _df1_id;
    j2c_array< double >  _df1_x;
    j2c_array< double >  _df1_y;
    int64_t __hpat_h5_dim_size_2_1;
    hsize_t* SSAValue0;
    double SSAValue1;
    int64_t parallel_ir_array_temp__3_4_1;
    int64_t parfor_index_1_3;
    int64_t parallel_ir_save_array_len_1_3;
    bool SSAValue2;
    j2c_array< bool >  parallel_ir_new_array_name_3_1;
    bool parallel_ir_array_temp__24_6_1;
    int64_t parfor_index_1_7;
    double parallel_ir_array_temp__13_9_1;
    int64_t parallel_ir_save_array_len_1_7;
    double parallel_ir_reduction_output_7;
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
    int64_t __hpat_loop_start_3;
    int64_t __hpat_loop_end_3;
    int64_t __hpat_loop_div_3;
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
    _df_id = j2c_array<int64_t>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
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
    ret_2 = H5Dread(dataset_id_2, H5T_NATIVE_LLONG, mem_dataspace_2, space_id_2, xfer_plist_2, _df_id.getData());
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
    _df_x = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_2);
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
    ret_3 = H5Dread(dataset_id_3, H5T_NATIVE_DOUBLE, mem_dataspace_3, space_id_3, xfer_plist_3, _df_x.getData());
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
    _df_y = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_3);
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
    ret_4 = H5Dread(dataset_id_4, H5T_NATIVE_DOUBLE, mem_dataspace_4, space_id_4, xfer_plist_4, _df_y.getData());
    assert(ret_4 != -1);
    ;
    ;
    H5Dclose(dataset_id_4);
    H5Fclose(file_id_4);
    ;
    parallel_ir_save_array_len_1_3 = __hpat_h5_dim_size_2_1;
    __hpat_dist_arr_div_4 = (parallel_ir_save_array_len_1_3) / (__hpat_num_pes);
    __hpat_dist_arr_start_4 = (__hpat_node_id) * (__hpat_dist_arr_div_4);
    __hpat_dist_arr_count_4 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_3-__hpat_node_id*__hpat_dist_arr_div_4 : __hpat_dist_arr_div_4);
    parallel_ir_new_array_name_3_1 = j2c_array<bool>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_4);
    __hpat_loop_div_3 = (parallel_ir_save_array_len_1_3) / (__hpat_num_pes);
    __hpat_loop_start_3 = ((__hpat_node_id) * (__hpat_loop_div_3)) + (1);
    __hpat_loop_end_3 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_3 : (__hpat_node_id+1)*__hpat_loop_div_3);
    for ( parfor_index_1_3 = __hpat_loop_start_3; parfor_index_1_3 <= (int64_t)__hpat_loop_end_3; parfor_index_1_3 += 1)
    {
        ;
        parallel_ir_array_temp__3_4_1 = _df_id.ARRAYELEM(((parfor_index_1_3) - (__hpat_loop_start_3)) + (1));
        SSAValue2 = (parallel_ir_array_temp__3_4_1) < (100);
        parallel_ir_array_temp__24_6_1 = SSAValue2;
        parallel_ir_new_array_name_3_1.ARRAYELEM(((parfor_index_1_3) - (__hpat_loop_start_3)) + (1)) = parallel_ir_array_temp__24_6_1;
    }
    ;
    _df_cond_e = parallel_ir_new_array_name_3_1;
    int _df_id_array_length_filter1 = _df_id.ARRAYLEN();
    int write_index_filter1 = 1;
    for (int index = 1 ; index < _df_id_array_length_filter1 + 1 ; index++)
    {
        if ( _df_cond_e.ARRAYELEM(index) )
        {
            _df_id.ARRAYELEM(write_index_filter1) =  _df_id.ARRAYELEM(index);
            _df_x.ARRAYELEM(write_index_filter1) =  _df_x.ARRAYELEM(index);
            _df_y.ARRAYELEM(write_index_filter1) =  _df_y.ARRAYELEM(index);
            write_index_filter1 = write_index_filter1 + 1;
        };
    };
    _df_id.dims[0] =  write_index_filter1 - 1;
    _df1_id = _df_id ;
    _df_x.dims[0] =  write_index_filter1 - 1;
    _df1_x = _df_x ;
    _df_y.dims[0] =  write_index_filter1 - 1;
    _df1_y = _df_y ;
    ;
    parallel_ir_save_array_len_1_7 = _df1_x.ARRAYSIZE(1);
    parallel_ir_reduction_output_7 = 0.0;
    for ( parfor_index_1_7 = 1; parfor_index_1_7 <= (int64_t)parallel_ir_save_array_len_1_7; parfor_index_1_7 += 1)
    {
        ;
        parallel_ir_array_temp__13_9_1 = _df1_x.ARRAYELEM(parfor_index_1_7);
        SSAValue3 = (parallel_ir_reduction_output_7) + (parallel_ir_array_temp__13_9_1);
        parallel_ir_reduction_output_7 = SSAValue3;
    }
    ;
    SSAValue1 = parallel_ir_reduction_output_7;
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
