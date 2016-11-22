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
} pppkernelscore_testp271;
typedef struct
{
    double f0;
    double f1;
    double f2;
} TupleFloat64Float64Float64;
typedef struct
{
} Basepvect;
typedef struct
{
    int64_t start;
    int64_t stop;
} UnitRangeInt64;
typedef struct
{
    int64_t f0;
    int64_t f1;
} TupleInt64Int64;
typedef struct
{
} ParallelAcceleratorAPIplog;
j2c_array< double >  _Base_vect(double X1, double X2, double X3);
double _ParallelAcceleratorAPI_log(double A);
j2c_array< double >  _Base_vect(double X1, double X2, double X3)
{
    TupleFloat64Float64Float64 X = {X1, X2, X3};
    Basepvect pselfp;
    int64_t ptempp;
    int64_t ppptempp_4p296;
    int64_t ppptempp_5p297;
    int64_t i;
    int64_t ppptempp_7p298;
    UnitRangeInt64 SSAValue0;
    int64_t SSAValue1;
    j2c_array< double >  SSAValue2;
    int64_t SSAValue3;
    TupleInt64Int64 SSAValue4;
    TupleInt64Int64 SSAValue5;
    TupleInt64Int64 SSAValue6;
    double SSAValue7;
    int64_t SSAValue8;
    int64_t SSAValue9;
    int64_t SSAValue10;
    int64_t SSAValue11;
    int64_t SSAValue12;
    int64_t SSAValue13;
    int64_t SSAValue14;
    SSAValue8 = 3;
    SSAValue10 = ((1) <= (SSAValue8)) ? (SSAValue8) : ((1) - (1));
    SSAValue1 = ((SSAValue10) - (1)) + (1);
    SSAValue2 = j2c_array<double>::new_j2c_array_1d(NULL, SSAValue1);
    ptempp = 1;
    ppptempp_4p296 = 1;
    ppptempp_5p297 = 0;
    label8 : ;
    if (!(!((ppptempp_5p297 == SSAValue1)))) goto label28;
    SSAValue3 = (ppptempp_5p297) + (1);
    ppptempp_5p297 = SSAValue3;
    SSAValue11 = ppptempp_4p296;
    SSAValue12 = (ppptempp_4p296) + (1);
    ppptempp_7p298 = 1;
    SSAValue13 = (1) + (1);
    i = SSAValue11;
    ppptempp_7p298 = SSAValue13;
    SSAValue14 = (2) + (1);
    ppptempp_4p296 = SSAValue12;
    ppptempp_7p298 = SSAValue14;
    SSAValue7 = ((double *)&X)[i - 1];
    SSAValue2.ARRAYELEM(ptempp) = SSAValue7;
    ptempp = (ptempp) + (1);
    label26 : ;
    goto label8;
    label28 : ;
    return SSAValue2;

}


double _ParallelAcceleratorAPI_log(double A)
{
    ParallelAcceleratorAPIplog pselfp;
    return log(A);;

}


void ppkernelscore_testp271(ASCIIString&  file_name, double * __restrict ret0)
{
    pppkernelscore_testp271 pselfp;
    j2c_array< int64_t >  _df_id;
    j2c_array< double >  _df_x;
    j2c_array< double >  _df_y;
    j2c_array< double >  points;
    int64_t N;
    double b;
    double exps;
    int64_t __hpat_h5_dim_size_1_1;
    int64_t SSAValue7;
    int64_t SSAValue9;
    hsize_t* SSAValue11;
    bool SSAValue17;
    int64_t SSAValue18;
    int64_t SSAValue20;
    int64_t SSAValue22;
    int64_t ppip273p279;
    double parallel_ir_reduction_input_1_1;
    double pptemp_neutral_valp280;
    int64_t i;
    j2c_array< double >  d;
    j2c_array< double >  SSAValue23;
    double SSAValue24;
    double SSAValue25;
    double SSAValue26;
    double SSAValue27;
    double SSAValue28;
    double SSAValue29;
    int64_t SSAValue30;
    int64_t SSAValue31;
    double SSAValue23pp1;
    int32_t SSAValue32;
    double SSAValue33;
    double SSAValue32pp2;
    double SSAValue34;
    double SSAValue35;
    double SSAValue38pp3;
    double SSAValue36;
    double parallel_ir_array_temp__6_12_1;
    int64_t parfor_index_1_11;
    int64_t parallel_ir_save_array_len_1_11;
    double SSAValue37;
    j2c_array< double >  parallel_ir_new_array_name_11_1;
    double parallel_ir_array_temp__23_14_1;
    double parallel_ir_array_temp_SSAValue24_16_1;
    int32_t SSAValue38;
    double SSAValue39;
    double parallel_ir_array_temp_SSAValue24_18_2;
    double parallel_ir_array_temp_SSAValue25_21_1;
    double SSAValue40;
    double parallel_ir_array_temp_SSAValue25_23_2;
    double parallel_ir_array_temp_SSAValue4_26_1;
    int64_t parfor_index_1_25;
    int64_t parallel_ir_save_array_len_1_25;
    double SSAValue41;
    double parallel_ir_array_temp_SSAValue4_28_2;
    double parallel_ir_array_temp__4_31_1;
    double parallel_ir_reduction_output_29;
    int64_t SSAValue42;
    int64_t SSAValue43;
    bool SSAValue44;
    bool SSAValue45;
    bool SSAValue46;
    bool SSAValue47;
    bool SSAValue48;
    bool SSAValue49;
    bool SSAValue50;
    bool SSAValue51;
    double SSAValue52;
    double SSAValue53;
    double SSAValue54;
    int64_t parfor_index_1_34;
    double parallel_ir_array_temp__4_36_1;
    int64_t parallel_ir_save_array_len_1_34;
    double parallel_ir_reduction_output_34;
    int64_t SSAValue0;
    int64_t SSAValue1;
    bool SSAValue2;
    bool SSAValue3;
    bool SSAValue4;
    bool SSAValue5;
    bool SSAValue6;
    bool SSAValue8;
    bool SSAValue10;
    bool SSAValue12;
    double SSAValue13;
    double SSAValue14;
    double SSAValue15;
    double parallel_ir_array_temp__4_39_1;
    int64_t parfor_index_1_38;
    int64_t parallel_ir_save_array_len_1_38;
    double SSAValue16;
    double parallel_ir_array_temp__4_41_2;
    double parallel_ir_array_temp_SSAValue39_43_1;
    double SSAValue19;
    double parallel_ir_array_temp_SSAValue39_45_2;
    double parallel_ir_array_temp_SSAValue40_49_1;
    double parallel_ir_reduction_output_47;
    double SSAValue21;
    int32_t __hpat_num_pes;
    int32_t __hpat_node_id;
    int64_t __hpat_dist_arr_start_1;
    int64_t __hpat_dist_arr_div_1;
    int64_t __hpat_dist_arr_count_1;
    int64_t __hpat_dist_arr_start_2;
    int64_t __hpat_dist_arr_div_2;
    int64_t __hpat_dist_arr_count_2;
    int64_t __hpat_loop_start_1;
    int64_t __hpat_loop_end_1;
    int64_t __hpat_loop_div_1;
    double __hpat_reduce_3;
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
    SSAValue11 = space_dims_1;;
    __hpat_h5_dim_size_1_1 = SSAValue11[1-1];
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
    _df_x = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_h5_dim_size_1_1);
    hsize_t CGen_HDF5_start_2[data_ndim_2];
    hsize_t CGen_HDF5_count_2[data_ndim_2];
    CGen_HDF5_start_2[0] = 0;
    CGen_HDF5_count_2[0] = space_dims_2[0];
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
    __hpat_dist_arr_div_2 = (__hpat_h5_dim_size_1_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_2 = (__hpat_node_id) * (__hpat_dist_arr_div_2);
    __hpat_dist_arr_count_2 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_1_1-__hpat_node_id*__hpat_dist_arr_div_2 : __hpat_dist_arr_div_2);
    _df_y = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_2);
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
    ret_3 = H5Dread(dataset_id_3, H5T_NATIVE_DOUBLE, mem_dataspace_3, space_id_3, xfer_plist_3, _df_y.getData());
    assert(ret_3 != -1);
    ;
    ;
    H5Dclose(dataset_id_3);
    H5Fclose(file_id_3);
    ;
    t1 = MPI_Wtime();
    points = _Base_vect(-1.0,2.0,5.0);
    N = points.ARRAYSIZE(1);
    b = 0.5;
    exps = 0.0;
    SSAValue7 = _df_x.ARRAYLEN();
    SSAValue17 = (1) <= (SSAValue7);
    SSAValue18 = (1) - (1);
    SSAValue9 = (SSAValue17) ? (SSAValue7) : (SSAValue18);
    SSAValue20 = (SSAValue9) - (1);
    SSAValue22 = (SSAValue20) + (1);
    __hpat_loop_div_1 = (SSAValue22) / (__hpat_num_pes);
    __hpat_loop_start_1 = ((__hpat_node_id) * (__hpat_loop_div_1)) + (1);
    __hpat_loop_end_1 = ((__hpat_node_id==__hpat_num_pes-1) ? SSAValue22 : (__hpat_node_id+1)*__hpat_loop_div_1);
    for ( ppip273p279 = __hpat_loop_start_1; ppip273p279 <= (int64_t)__hpat_loop_end_1; ppip273p279 += 1)
    {
        ;
        SSAValue30 = (ppip273p279) - (1);
        SSAValue31 = (SSAValue30) * (1);
        i = (1) + (SSAValue31);
        SSAValue23pp1 = _df_x.ARRAYELEM(i);
        parallel_ir_save_array_len_1_11 = points.ARRAYSIZE(1);
        parallel_ir_new_array_name_11_1 = j2c_array<double>::new_j2c_array_1d(NULL, parallel_ir_save_array_len_1_11);
        for ( parfor_index_1_11 = 1; parfor_index_1_11 <= (int64_t)parallel_ir_save_array_len_1_11; parfor_index_1_11 += 1)
        {
            ;
            parallel_ir_array_temp__6_12_1 = points.ARRAYELEM(parfor_index_1_11);
            SSAValue37 = (SSAValue23pp1) - (parallel_ir_array_temp__6_12_1);
            parallel_ir_array_temp__23_14_1 = SSAValue37;
            parallel_ir_array_temp_SSAValue24_16_1 = parallel_ir_array_temp__23_14_1;
            SSAValue38 = (int32_t)(2);
            pow(parallel_ir_array_temp_SSAValue24_16_1, SSAValue38);
            SSAValue39 = pow(parallel_ir_array_temp_SSAValue24_16_1, SSAValue38);
            parallel_ir_array_temp_SSAValue24_18_2 = SSAValue39;
            parallel_ir_array_temp_SSAValue25_21_1 = parallel_ir_array_temp_SSAValue24_18_2;
            SSAValue40 = -(parallel_ir_array_temp_SSAValue25_21_1);
            parallel_ir_array_temp_SSAValue25_23_2 = SSAValue40;
            parallel_ir_new_array_name_11_1.ARRAYELEM(parfor_index_1_11) = parallel_ir_array_temp_SSAValue25_23_2;
        }
        ;
        SSAValue23 = parallel_ir_new_array_name_11_1;
        SSAValue32 = (int32_t)(2);
        pow(b, SSAValue32);
        SSAValue25 = pow(b, SSAValue32);
        SSAValue33 = (double)2;
        SSAValue32pp2 = (SSAValue33) * (SSAValue25);
        parallel_ir_save_array_len_1_25 = SSAValue23.ARRAYSIZE(1);
        parallel_ir_reduction_output_29 = DBL_MAX;
        for ( parfor_index_1_25 = 1; parfor_index_1_25 <= (int64_t)parallel_ir_save_array_len_1_25; parfor_index_1_25 += 1)
        {
            ;
            parallel_ir_array_temp_SSAValue4_26_1 = SSAValue23.ARRAYELEM(parfor_index_1_25);
            SSAValue41 = (parallel_ir_array_temp_SSAValue4_26_1) / (SSAValue32pp2);
            parallel_ir_array_temp_SSAValue4_28_2 = SSAValue41;
            SSAValue23.ARRAYELEM(parfor_index_1_25) = parallel_ir_array_temp_SSAValue4_28_2;
            parallel_ir_array_temp__4_31_1 = parallel_ir_array_temp_SSAValue4_28_2;
            SSAValue42 = parallel_ir_array_temp__4_31_1;
            SSAValue43 = parallel_ir_reduction_output_29;
            SSAValue44 = (SSAValue43) < (0);
            SSAValue45 = (SSAValue42) < (0);
            SSAValue46 = !(SSAValue44);
            (SSAValue45) & (SSAValue46);
            SSAValue47 = (parallel_ir_array_temp__4_31_1) < (parallel_ir_reduction_output_29);
            SSAValue48 = (SSAValue45) & (SSAValue46);
            (SSAValue47) | (SSAValue48);
            SSAValue49 = (parallel_ir_array_temp__4_31_1) != (parallel_ir_array_temp__4_31_1);
            SSAValue50 = (parallel_ir_reduction_output_29) != (parallel_ir_reduction_output_29);
            SSAValue51 = (SSAValue47) | (SSAValue48);
            SSAValue52 = (SSAValue49) ? (parallel_ir_reduction_output_29) : (parallel_ir_array_temp__4_31_1);
            SSAValue53 = (SSAValue50) ? (parallel_ir_array_temp__4_31_1) : (parallel_ir_reduction_output_29);
            SSAValue54 = (SSAValue51) ? (SSAValue52) : (SSAValue53);
            parallel_ir_reduction_output_29 = SSAValue54;
        }
        ;
        d = SSAValue23;
        SSAValue24 = parallel_ir_reduction_output_29;
        SSAValue34 = (double)N;
        SSAValue35 = (b) * (SSAValue34);
        SSAValue26 = _ParallelAcceleratorAPI_log(SSAValue35);
        SSAValue27 = (SSAValue24) - (SSAValue26);
        parallel_ir_save_array_len_1_34 = d.ARRAYSIZE(1);
        parallel_ir_reduction_output_34 = DBL_MAX;
        for ( parfor_index_1_34 = 1; parfor_index_1_34 <= (int64_t)parallel_ir_save_array_len_1_34; parfor_index_1_34 += 1)
        {
            ;
            parallel_ir_array_temp__4_36_1 = d.ARRAYELEM(parfor_index_1_34);
            SSAValue0 = parallel_ir_array_temp__4_36_1;
            SSAValue1 = parallel_ir_reduction_output_34;
            SSAValue2 = (SSAValue1) < (0);
            SSAValue3 = (SSAValue0) < (0);
            SSAValue4 = !(SSAValue2);
            (SSAValue3) & (SSAValue4);
            SSAValue5 = (parallel_ir_array_temp__4_36_1) < (parallel_ir_reduction_output_34);
            SSAValue6 = (SSAValue3) & (SSAValue4);
            (SSAValue5) | (SSAValue6);
            SSAValue8 = (parallel_ir_array_temp__4_36_1) != (parallel_ir_array_temp__4_36_1);
            SSAValue10 = (parallel_ir_reduction_output_34) != (parallel_ir_reduction_output_34);
            SSAValue12 = (SSAValue5) | (SSAValue6);
            SSAValue13 = (SSAValue8) ? (parallel_ir_reduction_output_34) : (parallel_ir_array_temp__4_36_1);
            SSAValue14 = (SSAValue10) ? (parallel_ir_array_temp__4_36_1) : (parallel_ir_reduction_output_34);
            SSAValue15 = (SSAValue12) ? (SSAValue13) : (SSAValue14);
            parallel_ir_reduction_output_34 = SSAValue15;
        }
        ;
        SSAValue38pp3 = parallel_ir_reduction_output_34;
        parallel_ir_save_array_len_1_38 = d.ARRAYSIZE(1);
        parallel_ir_reduction_output_47 = 0.0;
        for ( parfor_index_1_38 = 1; parfor_index_1_38 <= (int64_t)parallel_ir_save_array_len_1_38; parfor_index_1_38 += 1)
        {
            ;
            parallel_ir_array_temp__4_39_1 = d.ARRAYELEM(parfor_index_1_38);
            SSAValue16 = (parallel_ir_array_temp__4_39_1) - (SSAValue38pp3);
            parallel_ir_array_temp__4_41_2 = SSAValue16;
            parallel_ir_array_temp_SSAValue39_43_1 = parallel_ir_array_temp__4_41_2;
            SSAValue19 = exp(parallel_ir_array_temp_SSAValue39_43_1);
            parallel_ir_array_temp_SSAValue39_45_2 = SSAValue19;
            parallel_ir_array_temp_SSAValue40_49_1 = parallel_ir_array_temp_SSAValue39_45_2;
            SSAValue21 = (parallel_ir_reduction_output_47) + (parallel_ir_array_temp_SSAValue40_49_1);
            parallel_ir_reduction_output_47 = SSAValue21;
        }
        ;
        SSAValue36 = parallel_ir_reduction_output_47;
        SSAValue28 = _ParallelAcceleratorAPI_log(SSAValue36);
        SSAValue29 = (SSAValue27) + (SSAValue28);
        exps = (exps) + (SSAValue29);
    }
    ;
    __hpat_reduce_3 = 0;
    MPI_Allreduce(&exps, &__hpat_reduce_3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
    exps = __hpat_reduce_3;
    t2 = MPI_Wtime();
    *ret0 = exps;
    return;

}


void ppkernelscore_testp271_unaliased(ASCIIString& __restrict file_name, double * __restrict ret0)
{
    pppkernelscore_testp271 pselfp;
    j2c_array< int64_t >  _df_id;
    j2c_array< double >  _df_x;
    j2c_array< double >  _df_y;
    j2c_array< double >  points;
    int64_t N;
    double b;
    double exps;
    int64_t __hpat_h5_dim_size_1_1;
    int64_t SSAValue7;
    int64_t SSAValue9;
    hsize_t* SSAValue11;
    bool SSAValue17;
    int64_t SSAValue18;
    int64_t SSAValue20;
    int64_t SSAValue22;
    int64_t ppip273p279;
    double parallel_ir_reduction_input_1_1;
    double pptemp_neutral_valp280;
    int64_t i;
    j2c_array< double >  d;
    j2c_array< double >  SSAValue23;
    double SSAValue24;
    double SSAValue25;
    double SSAValue26;
    double SSAValue27;
    double SSAValue28;
    double SSAValue29;
    int64_t SSAValue30;
    int64_t SSAValue31;
    double SSAValue23pp1;
    int32_t SSAValue32;
    double SSAValue33;
    double SSAValue32pp2;
    double SSAValue34;
    double SSAValue35;
    double SSAValue38pp3;
    double SSAValue36;
    double parallel_ir_array_temp__6_12_1;
    int64_t parfor_index_1_11;
    int64_t parallel_ir_save_array_len_1_11;
    double SSAValue37;
    j2c_array< double >  parallel_ir_new_array_name_11_1;
    double parallel_ir_array_temp__23_14_1;
    double parallel_ir_array_temp_SSAValue24_16_1;
    int32_t SSAValue38;
    double SSAValue39;
    double parallel_ir_array_temp_SSAValue24_18_2;
    double parallel_ir_array_temp_SSAValue25_21_1;
    double SSAValue40;
    double parallel_ir_array_temp_SSAValue25_23_2;
    double parallel_ir_array_temp_SSAValue4_26_1;
    int64_t parfor_index_1_25;
    int64_t parallel_ir_save_array_len_1_25;
    double SSAValue41;
    double parallel_ir_array_temp_SSAValue4_28_2;
    double parallel_ir_array_temp__4_31_1;
    double parallel_ir_reduction_output_29;
    int64_t SSAValue42;
    int64_t SSAValue43;
    bool SSAValue44;
    bool SSAValue45;
    bool SSAValue46;
    bool SSAValue47;
    bool SSAValue48;
    bool SSAValue49;
    bool SSAValue50;
    bool SSAValue51;
    double SSAValue52;
    double SSAValue53;
    double SSAValue54;
    int64_t parfor_index_1_34;
    double parallel_ir_array_temp__4_36_1;
    int64_t parallel_ir_save_array_len_1_34;
    double parallel_ir_reduction_output_34;
    int64_t SSAValue0;
    int64_t SSAValue1;
    bool SSAValue2;
    bool SSAValue3;
    bool SSAValue4;
    bool SSAValue5;
    bool SSAValue6;
    bool SSAValue8;
    bool SSAValue10;
    bool SSAValue12;
    double SSAValue13;
    double SSAValue14;
    double SSAValue15;
    double parallel_ir_array_temp__4_39_1;
    int64_t parfor_index_1_38;
    int64_t parallel_ir_save_array_len_1_38;
    double SSAValue16;
    double parallel_ir_array_temp__4_41_2;
    double parallel_ir_array_temp_SSAValue39_43_1;
    double SSAValue19;
    double parallel_ir_array_temp_SSAValue39_45_2;
    double parallel_ir_array_temp_SSAValue40_49_1;
    double parallel_ir_reduction_output_47;
    double SSAValue21;
    int32_t __hpat_num_pes;
    int32_t __hpat_node_id;
    int64_t __hpat_dist_arr_start_1;
    int64_t __hpat_dist_arr_div_1;
    int64_t __hpat_dist_arr_count_1;
    int64_t __hpat_dist_arr_start_2;
    int64_t __hpat_dist_arr_div_2;
    int64_t __hpat_dist_arr_count_2;
    int64_t __hpat_loop_start_1;
    int64_t __hpat_loop_end_1;
    int64_t __hpat_loop_div_1;
    double __hpat_reduce_3;
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
    SSAValue11 = space_dims_1;;
    __hpat_h5_dim_size_1_1 = SSAValue11[1-1];
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
    _df_x = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_h5_dim_size_1_1);
    hsize_t CGen_HDF5_start_2[data_ndim_2];
    hsize_t CGen_HDF5_count_2[data_ndim_2];
    CGen_HDF5_start_2[0] = 0;
    CGen_HDF5_count_2[0] = space_dims_2[0];
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
    __hpat_dist_arr_div_2 = (__hpat_h5_dim_size_1_1) / (__hpat_num_pes);
    __hpat_dist_arr_start_2 = (__hpat_node_id) * (__hpat_dist_arr_div_2);
    __hpat_dist_arr_count_2 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_1_1-__hpat_node_id*__hpat_dist_arr_div_2 : __hpat_dist_arr_div_2);
    _df_y = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_2);
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
    ret_3 = H5Dread(dataset_id_3, H5T_NATIVE_DOUBLE, mem_dataspace_3, space_id_3, xfer_plist_3, _df_y.getData());
    assert(ret_3 != -1);
    ;
    ;
    H5Dclose(dataset_id_3);
    H5Fclose(file_id_3);
    ;
    points = _Base_vect(-1.0,2.0,5.0);
    N = points.ARRAYSIZE(1);
    b = 0.5;
    exps = 0.0;
    SSAValue7 = _df_x.ARRAYLEN();
    SSAValue17 = (1) <= (SSAValue7);
    SSAValue18 = (1) - (1);
    SSAValue9 = (SSAValue17) ? (SSAValue7) : (SSAValue18);
    SSAValue20 = (SSAValue9) - (1);
    SSAValue22 = (SSAValue20) + (1);
    __hpat_loop_div_1 = (SSAValue22) / (__hpat_num_pes);
    __hpat_loop_start_1 = ((__hpat_node_id) * (__hpat_loop_div_1)) + (1);
    __hpat_loop_end_1 = ((__hpat_node_id==__hpat_num_pes-1) ? SSAValue22 : (__hpat_node_id+1)*__hpat_loop_div_1);
    for ( ppip273p279 = __hpat_loop_start_1; ppip273p279 <= (int64_t)__hpat_loop_end_1; ppip273p279 += 1)
    {
        ;
        SSAValue30 = (ppip273p279) - (1);
        SSAValue31 = (SSAValue30) * (1);
        i = (1) + (SSAValue31);
        SSAValue23pp1 = _df_x.ARRAYELEM(i);
        parallel_ir_save_array_len_1_11 = points.ARRAYSIZE(1);
        parallel_ir_new_array_name_11_1 = j2c_array<double>::new_j2c_array_1d(NULL, parallel_ir_save_array_len_1_11);
        for ( parfor_index_1_11 = 1; parfor_index_1_11 <= (int64_t)parallel_ir_save_array_len_1_11; parfor_index_1_11 += 1)
        {
            ;
            parallel_ir_array_temp__6_12_1 = points.ARRAYELEM(parfor_index_1_11);
            SSAValue37 = (SSAValue23pp1) - (parallel_ir_array_temp__6_12_1);
            parallel_ir_array_temp__23_14_1 = SSAValue37;
            parallel_ir_array_temp_SSAValue24_16_1 = parallel_ir_array_temp__23_14_1;
            SSAValue38 = (int32_t)(2);
            pow(parallel_ir_array_temp_SSAValue24_16_1, SSAValue38);
            SSAValue39 = pow(parallel_ir_array_temp_SSAValue24_16_1, SSAValue38);
            parallel_ir_array_temp_SSAValue24_18_2 = SSAValue39;
            parallel_ir_array_temp_SSAValue25_21_1 = parallel_ir_array_temp_SSAValue24_18_2;
            SSAValue40 = -(parallel_ir_array_temp_SSAValue25_21_1);
            parallel_ir_array_temp_SSAValue25_23_2 = SSAValue40;
            parallel_ir_new_array_name_11_1.ARRAYELEM(parfor_index_1_11) = parallel_ir_array_temp_SSAValue25_23_2;
        }
        ;
        SSAValue23 = parallel_ir_new_array_name_11_1;
        SSAValue32 = (int32_t)(2);
        pow(b, SSAValue32);
        SSAValue25 = pow(b, SSAValue32);
        SSAValue33 = (double)2;
        SSAValue32pp2 = (SSAValue33) * (SSAValue25);
        parallel_ir_save_array_len_1_25 = SSAValue23.ARRAYSIZE(1);
        parallel_ir_reduction_output_29 = DBL_MAX;
        for ( parfor_index_1_25 = 1; parfor_index_1_25 <= (int64_t)parallel_ir_save_array_len_1_25; parfor_index_1_25 += 1)
        {
            ;
            parallel_ir_array_temp_SSAValue4_26_1 = SSAValue23.ARRAYELEM(parfor_index_1_25);
            SSAValue41 = (parallel_ir_array_temp_SSAValue4_26_1) / (SSAValue32pp2);
            parallel_ir_array_temp_SSAValue4_28_2 = SSAValue41;
            SSAValue23.ARRAYELEM(parfor_index_1_25) = parallel_ir_array_temp_SSAValue4_28_2;
            parallel_ir_array_temp__4_31_1 = parallel_ir_array_temp_SSAValue4_28_2;
            SSAValue42 = parallel_ir_array_temp__4_31_1;
            SSAValue43 = parallel_ir_reduction_output_29;
            SSAValue44 = (SSAValue43) < (0);
            SSAValue45 = (SSAValue42) < (0);
            SSAValue46 = !(SSAValue44);
            (SSAValue45) & (SSAValue46);
            SSAValue47 = (parallel_ir_array_temp__4_31_1) < (parallel_ir_reduction_output_29);
            SSAValue48 = (SSAValue45) & (SSAValue46);
            (SSAValue47) | (SSAValue48);
            SSAValue49 = (parallel_ir_array_temp__4_31_1) != (parallel_ir_array_temp__4_31_1);
            SSAValue50 = (parallel_ir_reduction_output_29) != (parallel_ir_reduction_output_29);
            SSAValue51 = (SSAValue47) | (SSAValue48);
            SSAValue52 = (SSAValue49) ? (parallel_ir_reduction_output_29) : (parallel_ir_array_temp__4_31_1);
            SSAValue53 = (SSAValue50) ? (parallel_ir_array_temp__4_31_1) : (parallel_ir_reduction_output_29);
            SSAValue54 = (SSAValue51) ? (SSAValue52) : (SSAValue53);
            parallel_ir_reduction_output_29 = SSAValue54;
        }
        ;
        d = SSAValue23;
        SSAValue24 = parallel_ir_reduction_output_29;
        SSAValue34 = (double)N;
        SSAValue35 = (b) * (SSAValue34);
        SSAValue26 = _ParallelAcceleratorAPI_log(SSAValue35);
        SSAValue27 = (SSAValue24) - (SSAValue26);
        parallel_ir_save_array_len_1_34 = d.ARRAYSIZE(1);
        parallel_ir_reduction_output_34 = DBL_MAX;
        for ( parfor_index_1_34 = 1; parfor_index_1_34 <= (int64_t)parallel_ir_save_array_len_1_34; parfor_index_1_34 += 1)
        {
            ;
            parallel_ir_array_temp__4_36_1 = d.ARRAYELEM(parfor_index_1_34);
            SSAValue0 = parallel_ir_array_temp__4_36_1;
            SSAValue1 = parallel_ir_reduction_output_34;
            SSAValue2 = (SSAValue1) < (0);
            SSAValue3 = (SSAValue0) < (0);
            SSAValue4 = !(SSAValue2);
            (SSAValue3) & (SSAValue4);
            SSAValue5 = (parallel_ir_array_temp__4_36_1) < (parallel_ir_reduction_output_34);
            SSAValue6 = (SSAValue3) & (SSAValue4);
            (SSAValue5) | (SSAValue6);
            SSAValue8 = (parallel_ir_array_temp__4_36_1) != (parallel_ir_array_temp__4_36_1);
            SSAValue10 = (parallel_ir_reduction_output_34) != (parallel_ir_reduction_output_34);
            SSAValue12 = (SSAValue5) | (SSAValue6);
            SSAValue13 = (SSAValue8) ? (parallel_ir_reduction_output_34) : (parallel_ir_array_temp__4_36_1);
            SSAValue14 = (SSAValue10) ? (parallel_ir_array_temp__4_36_1) : (parallel_ir_reduction_output_34);
            SSAValue15 = (SSAValue12) ? (SSAValue13) : (SSAValue14);
            parallel_ir_reduction_output_34 = SSAValue15;
        }
        ;
        SSAValue38pp3 = parallel_ir_reduction_output_34;
        parallel_ir_save_array_len_1_38 = d.ARRAYSIZE(1);
        parallel_ir_reduction_output_47 = 0.0;
        for ( parfor_index_1_38 = 1; parfor_index_1_38 <= (int64_t)parallel_ir_save_array_len_1_38; parfor_index_1_38 += 1)
        {
            ;
            parallel_ir_array_temp__4_39_1 = d.ARRAYELEM(parfor_index_1_38);
            SSAValue16 = (parallel_ir_array_temp__4_39_1) - (SSAValue38pp3);
            parallel_ir_array_temp__4_41_2 = SSAValue16;
            parallel_ir_array_temp_SSAValue39_43_1 = parallel_ir_array_temp__4_41_2;
            SSAValue19 = exp(parallel_ir_array_temp_SSAValue39_43_1);
            parallel_ir_array_temp_SSAValue39_45_2 = SSAValue19;
            parallel_ir_array_temp_SSAValue40_49_1 = parallel_ir_array_temp_SSAValue39_45_2;
            SSAValue21 = (parallel_ir_reduction_output_47) + (parallel_ir_array_temp_SSAValue40_49_1);
            parallel_ir_reduction_output_47 = SSAValue21;
        }
        ;
        SSAValue36 = parallel_ir_reduction_output_47;
        SSAValue28 = _ParallelAcceleratorAPI_log(SSAValue36);
        SSAValue29 = (SSAValue27) + (SSAValue28);
        exps = (exps) + (SSAValue29);
    }
    ;
    __hpat_reduce_3 = 0;
    MPI_Allreduce(&exps, &__hpat_reduce_3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
    exps = __hpat_reduce_3;
    *ret0 = exps;
    return;

}


extern "C" void _ppkernelscore_testp271_unaliased_(int run_where, ASCIIString& __restrict file_name , double* __restrict ret0 , bool genMain = true)
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
        mainFile << "    _ppkernelscore_testp271_unaliased_(runwhere, file_name, &ret0, false);" << std::endl;
        mainFile << "    MPI_Finalize();" << std::endl;
        mainFile << "    return 0;" << std::endl;
        mainFile << "}" << std::endl;
        mainFile.close();
        std::ofstream mainFileSh(newMainSh.str());
        mainFileSh << "#!/bin/sh" << std::endl;
        mainFileSh << "mpiicpc -O3 -std=c++11  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << " -lhdf5  -lm " << std::endl;
        mainFileSh.close();
    }

    ppkernelscore_testp271_unaliased(file_name, ret0);
}


extern "C" void _ppkernelscore_testp271_(int run_where, ASCIIString&  file_name , double* __restrict ret0 , bool genMain = true)
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
        mainFile << "    _ppkernelscore_testp271_(runwhere, file_name, &ret0, false);" << std::endl;
        mainFile << "    MPI_Finalize();" << std::endl;
        mainFile << "    return 0;" << std::endl;
        mainFile << "}" << std::endl;
        mainFile.close();
        std::ofstream mainFileSh(newMainSh.str());
        mainFileSh << "#!/bin/sh" << std::endl;
        mainFileSh << "mpiicpc -O3 -std=c++11  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << " -lhdf5  -lm " << std::endl;
        mainFileSh.close();
    }

    ppkernelscore_testp271(file_name, ret0);
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
