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

#include <unordered_map>
unsigned main_count = 0;
typedef struct
{
} pppudf_testp271;
void ppudf_testp271(ASCIIString&  file_name, double * __restrict ret0)
{
    pppudf_testp271 pselfp;
    j2c_array< int64_t >  _df_id;
    j2c_array< double >  _df_x;
    j2c_array< double >  _df_y;
    j2c_array< double >  _df2_sx_e;
    j2c_array< double >  _df2_sy_e;
    j2c_array< int64_t >  _df2_id;
    j2c_array< double >  _df2_sx;
    j2c_array< double >  _df2_sy;
    int64_t __hpat_h5_dim_size_2_1;
    hsize_t* SSAValue4;
    double SSAValue9;
    double parallel_ir_array_temp__4_6_1;
    int64_t parfor_index_1_5;
    int64_t parallel_ir_save_array_len_1_5;
    double SSAValue10;
    double SSAValue11;
    double SSAValue0;
    double SSAValue1;
    j2c_array< double >  parallel_ir_new_array_name_5_1;
    double parallel_ir_array_temp__25_8_1;
    double parallel_ir_array_temp__5_10_1;
    double SSAValue2;
    double SSAValue3;
    double SSAValue5;
    double SSAValue6;
    j2c_array< double >  parallel_ir_new_array_name_9_1;
    double parallel_ir_array_temp__31_12_1;
    int64_t parfor_index_1_14;
    double parallel_ir_array_temp__12_16_1;
    int64_t parallel_ir_save_array_len_1_14;
    double parallel_ir_reduction_output_14;
    double SSAValue7;
    double SSAValue8;
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
    int64_t __hpat_loop_start_5;
    int64_t __hpat_loop_end_5;
    int64_t __hpat_loop_div_5;
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
    SSAValue4 = space_dims_2;;
    __hpat_h5_dim_size_2_1 = SSAValue4[1-1];
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
    parallel_ir_save_array_len_1_5 = __hpat_h5_dim_size_2_1;
    __hpat_dist_arr_div_4 = (parallel_ir_save_array_len_1_5) / (__hpat_num_pes);
    __hpat_dist_arr_start_4 = (__hpat_node_id) * (__hpat_dist_arr_div_4);
    __hpat_dist_arr_count_4 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_5-__hpat_node_id*__hpat_dist_arr_div_4 : __hpat_dist_arr_div_4);
    parallel_ir_new_array_name_5_1 = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_4);
    __hpat_dist_arr_div_5 = (parallel_ir_save_array_len_1_5) / (__hpat_num_pes);
    __hpat_dist_arr_start_5 = (__hpat_node_id) * (__hpat_dist_arr_div_5);
    __hpat_dist_arr_count_5 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_5-__hpat_node_id*__hpat_dist_arr_div_5 : __hpat_dist_arr_div_5);
    parallel_ir_new_array_name_9_1 = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_5);
    __hpat_loop_div_5 = (parallel_ir_save_array_len_1_5) / (__hpat_num_pes);
    __hpat_loop_start_5 = ((__hpat_node_id) * (__hpat_loop_div_5)) + (1);
    __hpat_loop_end_5 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_5 : (__hpat_node_id+1)*__hpat_loop_div_5);
    for ( parfor_index_1_5 = __hpat_loop_start_5; parfor_index_1_5 <= (int64_t)__hpat_loop_end_5; parfor_index_1_5 += 1)
    {
        ;
        parallel_ir_array_temp__4_6_1 = _df_x.ARRAYELEM(((parfor_index_1_5) - (__hpat_loop_start_5)) + (1));
        SSAValue10 = (double)2;
        SSAValue11 = SSAValue10;
        SSAValue0 = (SSAValue11) * (parallel_ir_array_temp__4_6_1);
        SSAValue1 = SSAValue0;
        parallel_ir_array_temp__25_8_1 = SSAValue1;
        parallel_ir_new_array_name_5_1.ARRAYELEM(((parfor_index_1_5) - (__hpat_loop_start_5)) + (1)) = parallel_ir_array_temp__25_8_1;
        parallel_ir_array_temp__5_10_1 = _df_y.ARRAYELEM(((parfor_index_1_5) - (__hpat_loop_start_5)) + (1));
        SSAValue2 = (double)2;
        SSAValue3 = SSAValue2;
        SSAValue5 = (SSAValue3) * (parallel_ir_array_temp__5_10_1);
        SSAValue6 = SSAValue5;
        parallel_ir_array_temp__31_12_1 = SSAValue6;
        parallel_ir_new_array_name_9_1.ARRAYELEM(((parfor_index_1_5) - (__hpat_loop_start_5)) + (1)) = parallel_ir_array_temp__31_12_1;
    }
    ;
    _df2_sx_e = parallel_ir_new_array_name_5_1;
    _df2_sy_e = parallel_ir_new_array_name_9_1;
    std::unordered_map<int,int> temp_map_1__df2_id;
    int *scount_1 = (int*)malloc(sizeof(int)*__hpat_num_pes);
    memset(scount_1, 0, sizeof(int)*__hpat_num_pes);
    int *s_ind_tmp_1 = (int*)malloc(sizeof(int)*__hpat_num_pes);
    memset(s_ind_tmp_1, 0, sizeof(int)*__hpat_num_pes);
    int rsize_1 = 0;
    int *rcount_1 = (int*)malloc(sizeof(int)*__hpat_num_pes);
    int *sdis_1 = (int*)malloc(sizeof(int)*__hpat_num_pes);
    int *rdis_1 = (int*)malloc(sizeof(int)*__hpat_num_pes);
    int agg_total_unique_keys_1 = 0;
    for (int i=1; i <= _df_id.ARRAYLEN(); i++)
    {
        if (temp_map_1__df2_id.find(_df_id.ARRAYELEM(i)) == temp_map_1__df2_id.end())
        {
            temp_map_1__df2_id[_df_id.ARRAYELEM(i)] = 1;
            int node_id = _df_id.ARRAYELEM(i) % __hpat_num_pes;
            scount_1[node_id]++;
            agg_total_unique_keys_1++;
        }
    }
    sdis_1[0]=0;
    for(int i=1;i < __hpat_num_pes;i++)
    {
        sdis_1[i] = scount_1[i-1] + sdis_1[i-1];
    }
    MPI_Alltoall(scount_1,1,MPI_INT,rcount_1,1,MPI_INT,MPI_COMM_WORLD);
    temp_map_1__df2_id.clear();
    rdis_1[0]=0;
    for(int i=1;i < __hpat_num_pes;i++)
    {
        rdis_1[i] = rcount_1[i-1] + rdis_1[i-1];
    }
    for(int i=0;i<__hpat_num_pes;i++)
    {
        rsize_1 = rsize_1 + rcount_1[i];
    }
    j2c_array< int64_t > _df_id_tmp_agg_1 = j2c_array< int64_t >::new_j2c_array_1d(NULL, agg_total_unique_keys_1);
    j2c_array< double > _df2_sx_e_tmp_agg_1 = j2c_array< double >::new_j2c_array_1d(NULL, agg_total_unique_keys_1);
    j2c_array< double > _df2_sy_e_tmp_agg_1 = j2c_array< double >::new_j2c_array_1d(NULL, agg_total_unique_keys_1);
    for(int i = 1; i <= _df_id.ARRAYLEN(); i++)
    {
        int64_t key = _df_id.ARRAYELEM(i);
        int node_id = key % __hpat_num_pes ;
        if (temp_map_1__df2_id.find(key) == temp_map_1__df2_id.end())
        {
            int agg_write_index_1 = sdis_1[node_id]+s_ind_tmp_1[node_id]+1 ;
            temp_map_1__df2_id[key] = agg_write_index_1;
            _df_id_tmp_agg_1.ARRAYELEM(agg_write_index_1) = key;
            _df2_sx_e_tmp_agg_1.ARRAYELEM(agg_write_index_1) = _df2_sx_e.ARRAYELEM(i);
            _df2_sy_e_tmp_agg_1.ARRAYELEM(agg_write_index_1) = _df2_sy_e.ARRAYELEM(i);
            s_ind_tmp_1[node_id]++;
        }
        else
        {
            int current_write_index1 = temp_map_1__df2_id[key];
            _df2_sx_e_tmp_agg_1.ARRAYELEM(current_write_index1) += _df2_sx_e.ARRAYELEM(i);
            _df2_sy_e_tmp_agg_1.ARRAYELEM(current_write_index1) += _df2_sy_e.ARRAYELEM(i);
        }
    }
    temp_map_1__df2_id.clear();
    j2c_array< int64_t > rbuf_1__df_id = j2c_array< int64_t >::new_j2c_array_1d(NULL, rsize_1);
    MPI_Alltoallv(_df_id_tmp_agg_1.getData(), scount_1, sdis_1, MPI_INT64_T,
        rbuf_1__df_id.getData(), rcount_1, rdis_1, MPI_INT64_T, MPI_COMM_WORLD);
    j2c_array< double > rbuf_1__df2_sx_e = j2c_array< double >::new_j2c_array_1d(NULL, rsize_1);
    MPI_Alltoallv(_df2_sx_e_tmp_agg_1.getData(), scount_1, sdis_1, MPI_DOUBLE,
        rbuf_1__df2_sx_e.getData(), rcount_1, rdis_1, MPI_DOUBLE, MPI_COMM_WORLD);
    j2c_array< double > rbuf_1__df2_sy_e = j2c_array< double >::new_j2c_array_1d(NULL, rsize_1);
    MPI_Alltoallv(_df2_sy_e_tmp_agg_1.getData(), scount_1, sdis_1, MPI_DOUBLE,
        rbuf_1__df2_sy_e.getData(), rcount_1, rdis_1, MPI_DOUBLE, MPI_COMM_WORLD);
    _df2_id = j2c_array< int64_t >::new_j2c_array_1d(NULL, rbuf_1__df_id.ARRAYLEN());
    _df2_sx = j2c_array< double >::new_j2c_array_1d(NULL, rbuf_1__df_id.ARRAYLEN());
    _df2_sy = j2c_array< double >::new_j2c_array_1d(NULL, rbuf_1__df_id.ARRAYLEN());
    int agg_write_index_1 = 1;
    for(int i = 1 ; i < rbuf_1__df_id.ARRAYLEN() + 1 ; i++)
    {
        int64_t key = rbuf_1__df_id.ARRAYELEM(i);
        if (temp_map_1__df2_id.find(key) == temp_map_1__df2_id.end())
        {
            temp_map_1__df2_id[key] = agg_write_index_1 ;
            _df2_id.ARRAYELEM(agg_write_index_1) = key;
            _df2_sx.ARRAYELEM(agg_write_index_1) = rbuf_1__df2_sx_e.ARRAYELEM(i);
            _df2_sy.ARRAYELEM(agg_write_index_1) = rbuf_1__df2_sy_e.ARRAYELEM(i);
            agg_write_index_1++;
        }
        else
        {
            int current_write_index1 = temp_map_1__df2_id[key];
            _df2_sx.ARRAYELEM(current_write_index1) +=  rbuf_1__df2_sx_e.ARRAYELEM(i);
            _df2_sy.ARRAYELEM(current_write_index1) +=  rbuf_1__df2_sy_e.ARRAYELEM(i);
        }
    }
    int counter_agg_1 = temp_map_1__df2_id.size();
    _df2_id.dims[0] = counter_agg_1;
    _df2_sx.dims[0] = counter_agg_1;
    _df2_sy.dims[0] = counter_agg_1;
    ;
    parallel_ir_save_array_len_1_14 = _df2_sx.ARRAYSIZE(1);
    parallel_ir_reduction_output_14 = 0.0;
    for ( parfor_index_1_14 = 1; parfor_index_1_14 <= (int64_t)parallel_ir_save_array_len_1_14; parfor_index_1_14 += 1)
    {
        ;
        parallel_ir_array_temp__12_16_1 = _df2_sx.ARRAYELEM(parfor_index_1_14);
        SSAValue7 = (parallel_ir_reduction_output_14) + (parallel_ir_array_temp__12_16_1);
        SSAValue8 = SSAValue7;
        parallel_ir_reduction_output_14 = SSAValue8;
    }
    ;
    SSAValue9 = parallel_ir_reduction_output_14;
    *ret0 = SSAValue9;
    return;

}


void ppudf_testp271_unaliased(ASCIIString& __restrict file_name, double * __restrict ret0)
{
    pppudf_testp271 pselfp;
    j2c_array< int64_t >  _df_id;
    j2c_array< double >  _df_x;
    j2c_array< double >  _df_y;
    j2c_array< double >  _df2_sx_e;
    j2c_array< double >  _df2_sy_e;
    j2c_array< int64_t >  _df2_id;
    j2c_array< double >  _df2_sx;
    j2c_array< double >  _df2_sy;
    int64_t __hpat_h5_dim_size_2_1;
    hsize_t* SSAValue4;
    double SSAValue9;
    double parallel_ir_array_temp__4_6_1;
    int64_t parfor_index_1_5;
    int64_t parallel_ir_save_array_len_1_5;
    double SSAValue10;
    double SSAValue11;
    double SSAValue0;
    double SSAValue1;
    j2c_array< double >  parallel_ir_new_array_name_5_1;
    double parallel_ir_array_temp__25_8_1;
    double parallel_ir_array_temp__5_10_1;
    double SSAValue2;
    double SSAValue3;
    double SSAValue5;
    double SSAValue6;
    j2c_array< double >  parallel_ir_new_array_name_9_1;
    double parallel_ir_array_temp__31_12_1;
    int64_t parfor_index_1_14;
    double parallel_ir_array_temp__12_16_1;
    int64_t parallel_ir_save_array_len_1_14;
    double parallel_ir_reduction_output_14;
    double SSAValue7;
    double SSAValue8;
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
    int64_t __hpat_loop_start_5;
    int64_t __hpat_loop_end_5;
    int64_t __hpat_loop_div_5;
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
    SSAValue4 = space_dims_2;;
    __hpat_h5_dim_size_2_1 = SSAValue4[1-1];
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
    parallel_ir_save_array_len_1_5 = __hpat_h5_dim_size_2_1;
    __hpat_dist_arr_div_4 = (parallel_ir_save_array_len_1_5) / (__hpat_num_pes);
    __hpat_dist_arr_start_4 = (__hpat_node_id) * (__hpat_dist_arr_div_4);
    __hpat_dist_arr_count_4 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_5-__hpat_node_id*__hpat_dist_arr_div_4 : __hpat_dist_arr_div_4);
    parallel_ir_new_array_name_5_1 = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_4);
    __hpat_dist_arr_div_5 = (parallel_ir_save_array_len_1_5) / (__hpat_num_pes);
    __hpat_dist_arr_start_5 = (__hpat_node_id) * (__hpat_dist_arr_div_5);
    __hpat_dist_arr_count_5 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_5-__hpat_node_id*__hpat_dist_arr_div_5 : __hpat_dist_arr_div_5);
    parallel_ir_new_array_name_9_1 = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_5);
    __hpat_loop_div_5 = (parallel_ir_save_array_len_1_5) / (__hpat_num_pes);
    __hpat_loop_start_5 = ((__hpat_node_id) * (__hpat_loop_div_5)) + (1);
    __hpat_loop_end_5 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_5 : (__hpat_node_id+1)*__hpat_loop_div_5);
    for ( parfor_index_1_5 = __hpat_loop_start_5; parfor_index_1_5 <= (int64_t)__hpat_loop_end_5; parfor_index_1_5 += 1)
    {
        ;
        parallel_ir_array_temp__4_6_1 = _df_x.ARRAYELEM(((parfor_index_1_5) - (__hpat_loop_start_5)) + (1));
        SSAValue10 = (double)2;
        SSAValue11 = SSAValue10;
        SSAValue0 = (SSAValue11) * (parallel_ir_array_temp__4_6_1);
        SSAValue1 = SSAValue0;
        parallel_ir_array_temp__25_8_1 = SSAValue1;
        parallel_ir_new_array_name_5_1.ARRAYELEM(((parfor_index_1_5) - (__hpat_loop_start_5)) + (1)) = parallel_ir_array_temp__25_8_1;
        parallel_ir_array_temp__5_10_1 = _df_y.ARRAYELEM(((parfor_index_1_5) - (__hpat_loop_start_5)) + (1));
        SSAValue2 = (double)2;
        SSAValue3 = SSAValue2;
        SSAValue5 = (SSAValue3) * (parallel_ir_array_temp__5_10_1);
        SSAValue6 = SSAValue5;
        parallel_ir_array_temp__31_12_1 = SSAValue6;
        parallel_ir_new_array_name_9_1.ARRAYELEM(((parfor_index_1_5) - (__hpat_loop_start_5)) + (1)) = parallel_ir_array_temp__31_12_1;
    }
    ;
    _df2_sx_e = parallel_ir_new_array_name_5_1;
    _df2_sy_e = parallel_ir_new_array_name_9_1;
    std::unordered_map<int,int> temp_map_1__df2_id;
    int *scount_1 = (int*)malloc(sizeof(int)*__hpat_num_pes);
    memset(scount_1, 0, sizeof(int)*__hpat_num_pes);
    int *s_ind_tmp_1 = (int*)malloc(sizeof(int)*__hpat_num_pes);
    memset(s_ind_tmp_1, 0, sizeof(int)*__hpat_num_pes);
    int rsize_1 = 0;
    int *rcount_1 = (int*)malloc(sizeof(int)*__hpat_num_pes);
    int *sdis_1 = (int*)malloc(sizeof(int)*__hpat_num_pes);
    int *rdis_1 = (int*)malloc(sizeof(int)*__hpat_num_pes);
    int agg_total_unique_keys_1 = 0;
    for (int i=1; i <= _df_id.ARRAYLEN(); i++)
    {
        if (temp_map_1__df2_id.find(_df_id.ARRAYELEM(i)) == temp_map_1__df2_id.end())
        {
            temp_map_1__df2_id[_df_id.ARRAYELEM(i)] = 1;
            int node_id = _df_id.ARRAYELEM(i) % __hpat_num_pes;
            scount_1[node_id]++;
            agg_total_unique_keys_1++;
        }
    }
    sdis_1[0]=0;
    for(int i=1;i < __hpat_num_pes;i++)
    {
        sdis_1[i] = scount_1[i-1] + sdis_1[i-1];
    }
    MPI_Alltoall(scount_1,1,MPI_INT,rcount_1,1,MPI_INT,MPI_COMM_WORLD);
    temp_map_1__df2_id.clear();
    rdis_1[0]=0;
    for(int i=1;i < __hpat_num_pes;i++)
    {
        rdis_1[i] = rcount_1[i-1] + rdis_1[i-1];
    }
    for(int i=0;i<__hpat_num_pes;i++)
    {
        rsize_1 = rsize_1 + rcount_1[i];
    }
    j2c_array< int64_t > _df_id_tmp_agg_1 = j2c_array< int64_t >::new_j2c_array_1d(NULL, agg_total_unique_keys_1);
    j2c_array< double > _df2_sx_e_tmp_agg_1 = j2c_array< double >::new_j2c_array_1d(NULL, agg_total_unique_keys_1);
    j2c_array< double > _df2_sy_e_tmp_agg_1 = j2c_array< double >::new_j2c_array_1d(NULL, agg_total_unique_keys_1);
    for(int i = 1; i <= _df_id.ARRAYLEN(); i++)
    {
        int64_t key = _df_id.ARRAYELEM(i);
        int node_id = key % __hpat_num_pes ;
        if (temp_map_1__df2_id.find(key) == temp_map_1__df2_id.end())
        {
            int agg_write_index_1 = sdis_1[node_id]+s_ind_tmp_1[node_id]+1 ;
            temp_map_1__df2_id[key] = agg_write_index_1;
            _df_id_tmp_agg_1.ARRAYELEM(agg_write_index_1) = key;
            _df2_sx_e_tmp_agg_1.ARRAYELEM(agg_write_index_1) = _df2_sx_e.ARRAYELEM(i);
            _df2_sy_e_tmp_agg_1.ARRAYELEM(agg_write_index_1) = _df2_sy_e.ARRAYELEM(i);
            s_ind_tmp_1[node_id]++;
        }
        else
        {
            int current_write_index1 = temp_map_1__df2_id[key];
            _df2_sx_e_tmp_agg_1.ARRAYELEM(current_write_index1) += _df2_sx_e.ARRAYELEM(i);
            _df2_sy_e_tmp_agg_1.ARRAYELEM(current_write_index1) += _df2_sy_e.ARRAYELEM(i);
        }
    }
    temp_map_1__df2_id.clear();
    j2c_array< int64_t > rbuf_1__df_id = j2c_array< int64_t >::new_j2c_array_1d(NULL, rsize_1);
    MPI_Alltoallv(_df_id_tmp_agg_1.getData(), scount_1, sdis_1, MPI_INT64_T,
        rbuf_1__df_id.getData(), rcount_1, rdis_1, MPI_INT64_T, MPI_COMM_WORLD);
    j2c_array< double > rbuf_1__df2_sx_e = j2c_array< double >::new_j2c_array_1d(NULL, rsize_1);
    MPI_Alltoallv(_df2_sx_e_tmp_agg_1.getData(), scount_1, sdis_1, MPI_DOUBLE,
        rbuf_1__df2_sx_e.getData(), rcount_1, rdis_1, MPI_DOUBLE, MPI_COMM_WORLD);
    j2c_array< double > rbuf_1__df2_sy_e = j2c_array< double >::new_j2c_array_1d(NULL, rsize_1);
    MPI_Alltoallv(_df2_sy_e_tmp_agg_1.getData(), scount_1, sdis_1, MPI_DOUBLE,
        rbuf_1__df2_sy_e.getData(), rcount_1, rdis_1, MPI_DOUBLE, MPI_COMM_WORLD);
    _df2_id = j2c_array< int64_t >::new_j2c_array_1d(NULL, rbuf_1__df_id.ARRAYLEN());
    _df2_sx = j2c_array< double >::new_j2c_array_1d(NULL, rbuf_1__df_id.ARRAYLEN());
    _df2_sy = j2c_array< double >::new_j2c_array_1d(NULL, rbuf_1__df_id.ARRAYLEN());
    int agg_write_index_1 = 1;
    for(int i = 1 ; i < rbuf_1__df_id.ARRAYLEN() + 1 ; i++)
    {
        int64_t key = rbuf_1__df_id.ARRAYELEM(i);
        if (temp_map_1__df2_id.find(key) == temp_map_1__df2_id.end())
        {
            temp_map_1__df2_id[key] = agg_write_index_1 ;
            _df2_id.ARRAYELEM(agg_write_index_1) = key;
            _df2_sx.ARRAYELEM(agg_write_index_1) = rbuf_1__df2_sx_e.ARRAYELEM(i);
            _df2_sy.ARRAYELEM(agg_write_index_1) = rbuf_1__df2_sy_e.ARRAYELEM(i);
            agg_write_index_1++;
        }
        else
        {
            int current_write_index1 = temp_map_1__df2_id[key];
            _df2_sx.ARRAYELEM(current_write_index1) +=  rbuf_1__df2_sx_e.ARRAYELEM(i);
            _df2_sy.ARRAYELEM(current_write_index1) +=  rbuf_1__df2_sy_e.ARRAYELEM(i);
        }
    }
    int counter_agg_1 = temp_map_1__df2_id.size();
    _df2_id.dims[0] = counter_agg_1;
    _df2_sx.dims[0] = counter_agg_1;
    _df2_sy.dims[0] = counter_agg_1;
    ;
    parallel_ir_save_array_len_1_14 = _df2_sx.ARRAYSIZE(1);
    parallel_ir_reduction_output_14 = 0.0;
    for ( parfor_index_1_14 = 1; parfor_index_1_14 <= (int64_t)parallel_ir_save_array_len_1_14; parfor_index_1_14 += 1)
    {
        ;
        parallel_ir_array_temp__12_16_1 = _df2_sx.ARRAYELEM(parfor_index_1_14);
        SSAValue7 = (parallel_ir_reduction_output_14) + (parallel_ir_array_temp__12_16_1);
        SSAValue8 = SSAValue7;
        parallel_ir_reduction_output_14 = SSAValue8;
    }
    ;
    SSAValue9 = parallel_ir_reduction_output_14;
    *ret0 = SSAValue9;
    return;

}


extern "C" void _ppudf_testp271_unaliased_(int run_where, ASCIIString& __restrict file_name , double* __restrict ret0 , bool genMain = true)
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
        mainFile << "    _ppudf_testp271_unaliased_(runwhere, file_name, &ret0, false);" << std::endl;
        mainFile << "    MPI_Finalize();" << std::endl;
        mainFile << "    return 0;" << std::endl;
        mainFile << "}" << std::endl;
        mainFile.close();
        std::ofstream mainFileSh(newMainSh.str());
        mainFileSh << "#!/bin/sh" << std::endl;
        mainFileSh << "mpiicpc -O3 -std=c++11  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << " -lhdf5  -lm " << std::endl;
        mainFileSh.close();
    }

    ppudf_testp271_unaliased(file_name, ret0);
}


extern "C" void _ppudf_testp271_(int run_where, ASCIIString&  file_name , double* __restrict ret0 , bool genMain = true)
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
        mainFile << "    _ppudf_testp271_(runwhere, file_name, &ret0, false);" << std::endl;
        mainFile << "    MPI_Finalize();" << std::endl;
        mainFile << "    return 0;" << std::endl;
        mainFile << "}" << std::endl;
        mainFile.close();
        std::ofstream mainFileSh(newMainSh.str());
        mainFileSh << "#!/bin/sh" << std::endl;
        mainFileSh << "mpiicpc -O3 -std=c++11  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << " -lhdf5  -lm " << std::endl;
        mainFileSh.close();
    }

    ppudf_testp271(file_name, ret0);
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
