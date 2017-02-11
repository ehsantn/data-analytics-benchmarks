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
} ppplogistic_regressionp271;
typedef struct
{
	int64_t f0;
	int64_t f1;
} TupleInt64Int64;
void pplogistic_regressionp271(int64_t iterations, ASCIIString&  file_name,  j2c_array< double >  * __restrict ret0)
{
	ppplogistic_regressionp271 pselfp;
	j2c_array< double >  points;
	j2c_array< double >  responses;
	j2c_array< double >  labels;
	j2c_array< double >  w;
	int64_t ppptempp_11p280;
	int64_t __hpat_h5_dim_size_1_2;
	int64_t __hpat_h5_dim_size_1_1;
	int64_t __hpat_h5_dim_size_2_1;
	double SSAValue40pp2;
	j2c_array< double >  SSAValue12;
	j2c_array< double >  SSAValue18;
	int64_t SSAValue0;
	hsize_t* SSAValue2;
	hsize_t* SSAValue3;
	TupleInt64Int64 SSAValue4;
	j2c_array< double >  SSAValue5;
	bool SSAValue6;
	int64_t SSAValue7;
	int64_t SSAValue8;
	bool SSAValue9;
	bool SSAValue10;
	j2c_array< double >  SSAValue11;
	j2c_array< double >  SSAValue13;
	int64_t parfor_index_1_53;
	int64_t parfor_index_2_53;
	double SSAValue14;
	double parallel_ir_array_temp__9_62_1;
	int64_t parfor_index_1_61;
	int64_t parfor_index_2_61;
	double SSAValue15;
	j2c_array< double >  parallel_ir_new_array_name_61_1;
	double parallel_ir_array_temp_SSAValue53_67_1;
	double SSAValue16;
	double SSAValue17;
	double SSAValue19;
	double SSAValue20;
	double SSAValue21;
	double parallel_ir_array_temp__9_88_1;
	double SSAValue22;
	double parallel_ir_array_temp__10_92_1;
	double parallel_ir_array_temp_SSAValue61_93_1;
	int64_t parfor_index_1_91;
	int64_t parfor_index_2_91;
	double SSAValue1;
	j2c_array< double >  parallel_ir_new_array_name_91_1;
	int64_t _dist_parfor_94_index1;
	int64_t _dist_parfor_95_index2;
	int64_t _dist_loop_96_index;
	double _dist_gemm_tmp1_97;
	double _dist_gemm_tmp2_98;
	double _dist_gemm_tmp3_99;
	int64_t _dist_parfor_101_index1;
	int64_t _dist_parfor_102_index2;
	int64_t _dist_loop_103_index;
	double _dist_gemm_tmp1_104;
	double _dist_gemm_tmp2_105;
	double _dist_gemm_tmp3_106;
	int32_t __hpat_num_pes;
	int32_t __hpat_node_id;
	int64_t __hpat_dist_arr_start_108;
	int64_t __hpat_dist_arr_div_108;
	int64_t __hpat_dist_arr_count_108;
	int64_t __hpat_dist_arr_start_109;
	int64_t __hpat_dist_arr_div_109;
	int64_t __hpat_dist_arr_count_109;
	int64_t __hpat_dist_arr_start_110;
	int64_t __hpat_dist_arr_div_110;
	int64_t __hpat_dist_arr_count_110;
	TupleInt64Int64 __hpat_dist_tup_var_110;
	int64_t __hpat_loop_start_100;
	int64_t __hpat_loop_end_100;
	int64_t __hpat_loop_div_100;
	j2c_array< double >  __hpat_reduce_111;
	;;
	MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
	MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
	SSAValue40pp2 = (double)1;
	SSAValue6 = (1) <= (iterations);
	SSAValue7 = (1) - (1);
	SSAValue0 = (SSAValue6) ? (iterations) : (SSAValue7);
	SSAValue8 = (SSAValue0) + (1);
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
	dataset_id_1 = H5Dopen2(file_id_1, "/points", H5P_DEFAULT);
	assert(dataset_id_1 != -1);
	;
	hid_t space_id_1 = H5Dget_space(dataset_id_1);
	assert(space_id_1 != -1);
	hsize_t data_ndim_1 = H5Sget_simple_extent_ndims(space_id_1);
	hsize_t space_dims_1[data_ndim_1];
	H5Sget_simple_extent_dims(space_id_1, space_dims_1, NULL);
	SSAValue2 = space_dims_1;;
	__hpat_h5_dim_size_1_2 = SSAValue2[2-1];
	__hpat_h5_dim_size_1_1 = SSAValue2[1-1];
	__hpat_dist_arr_div_108 = (__hpat_h5_dim_size_1_1) / (__hpat_num_pes);
	__hpat_dist_arr_start_108 = (__hpat_node_id) * (__hpat_dist_arr_div_108);
	__hpat_dist_arr_count_108 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_1_1-__hpat_node_id*__hpat_dist_arr_div_108 : __hpat_dist_arr_div_108);
	points = j2c_array<double>::new_j2c_array_2d(NULL, __hpat_h5_dim_size_1_2, __hpat_dist_arr_count_108);
	hsize_t CGen_HDF5_start_1[data_ndim_1];
	hsize_t CGen_HDF5_count_1[data_ndim_1];
	CGen_HDF5_start_1[0] = __hpat_dist_arr_start_108;
	CGen_HDF5_count_1[0] = __hpat_dist_arr_count_108;
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
	ret_1 = H5Dread(dataset_id_1, H5T_NATIVE_DOUBLE, mem_dataspace_1, space_id_1, xfer_plist_1, points.getData());
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
	dataset_id_2 = H5Dopen2(file_id_2, "/responses", H5P_DEFAULT);
	assert(dataset_id_2 != -1);
	;
	hid_t space_id_2 = H5Dget_space(dataset_id_2);
	assert(space_id_2 != -1);
	hsize_t data_ndim_2 = H5Sget_simple_extent_ndims(space_id_2);
	hsize_t space_dims_2[data_ndim_2];
	H5Sget_simple_extent_dims(space_id_2, space_dims_2, NULL);
	SSAValue3 = space_dims_2;;
	__hpat_h5_dim_size_2_1 = SSAValue3[1-1];
	__hpat_dist_arr_div_109 = (__hpat_h5_dim_size_2_1) / (__hpat_num_pes);
	__hpat_dist_arr_start_109 = (__hpat_node_id) * (__hpat_dist_arr_div_109);
	__hpat_dist_arr_count_109 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_2_1-__hpat_node_id*__hpat_dist_arr_div_109 : __hpat_dist_arr_div_109);
	responses = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_109);
	hsize_t CGen_HDF5_start_2[data_ndim_2];
	hsize_t CGen_HDF5_count_2[data_ndim_2];
	CGen_HDF5_start_2[0] = __hpat_dist_arr_start_109;
	CGen_HDF5_count_2[0] = __hpat_dist_arr_count_109;
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
	ret_2 = H5Dread(dataset_id_2, H5T_NATIVE_DOUBLE, mem_dataspace_2, space_id_2, xfer_plist_2, responses.getData());
	assert(ret_2 != -1);
	;
	;
	H5Dclose(dataset_id_2);
	H5Fclose(file_id_2);
	;
	SSAValue4 = TupleInt64Int64{1, __hpat_h5_dim_size_1_1};
	__hpat_dist_arr_div_110 = (__hpat_h5_dim_size_1_1) / (__hpat_num_pes);
	__hpat_dist_arr_start_110 = (__hpat_node_id) * (__hpat_dist_arr_div_110);
	__hpat_dist_arr_count_110 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_1_1-__hpat_node_id*__hpat_dist_arr_div_110 : __hpat_dist_arr_div_110);
	__hpat_dist_tup_var_110 = TupleInt64Int64{1, __hpat_dist_arr_count_110};
	labels = responses.reshape(__hpat_dist_tup_var_110.f0,__hpat_dist_tup_var_110.f1);
	;
	0; double __hpat_t1 = MPI_Wtime();
	ppptempp_11p280 = 1;
	SSAValue5 = j2c_array<double>::new_j2c_array_2d(NULL, 1, __hpat_h5_dim_size_1_2);
	for ( parfor_index_2_53 = 1; parfor_index_2_53 <= (int64_t)__hpat_h5_dim_size_1_2; parfor_index_2_53 += 1)
	{
		for ( parfor_index_1_53 = 1; parfor_index_1_53 <= (int64_t)1; parfor_index_1_53 += 1)
		{
			;
			SSAValue14 = (SSAValue40pp2) - (0.5);
			SSAValue5.ARRAYELEM(parfor_index_1_53,parfor_index_2_53) = SSAValue14;
		}
	}
	;
	w = SSAValue5;
	while (1)
	{
		;
		SSAValue9 = (ppptempp_11p280 == SSAValue8);
		SSAValue10 = !(SSAValue9);
		if (!(SSAValue10)) break;
		ppptempp_11p280 = (ppptempp_11p280) + (1);
		__hpat_loop_div_100 = (__hpat_h5_dim_size_1_1) / (__hpat_num_pes);
		__hpat_loop_start_100 = ((__hpat_node_id) * (__hpat_loop_div_100)) + (1);
		__hpat_loop_end_100 = ((__hpat_node_id==__hpat_num_pes-1) ? __hpat_h5_dim_size_1_1 : (__hpat_node_id+1)*__hpat_loop_div_100);
		SSAValue13 = j2c_array<double>::new_j2c_array_2d(NULL, 1, __hpat_h5_dim_size_1_2);
		memset(SSAValue13.data, 0, sizeof(double)*(__hpat_h5_dim_size_1_2) * (1));
		#pragma simd
		for ( _dist_parfor_94_index1 = __hpat_loop_start_100; _dist_parfor_94_index1 <= (int64_t)__hpat_loop_end_100; _dist_parfor_94_index1 += 1)
		{
			for ( _dist_parfor_95_index2 = 1; _dist_parfor_95_index2 <= (int64_t)1; _dist_parfor_95_index2 += 1)
			{
				;
				_dist_gemm_tmp3_99 = 0;
				for(int64_t _dist_loop_96_index = 1; _dist_loop_96_index <= __hpat_h5_dim_size_1_2; _dist_loop_96_index += 1)
				{
					;
					_dist_gemm_tmp1_97 = w.ARRAYELEM(_dist_parfor_95_index2,_dist_loop_96_index);
					_dist_gemm_tmp2_98 = points.ARRAYELEM(_dist_loop_96_index,((_dist_parfor_94_index1) - (__hpat_loop_start_100)) + (1));
					_dist_gemm_tmp3_99 = (_dist_gemm_tmp3_99) + ((_dist_gemm_tmp1_97) * (_dist_gemm_tmp2_98));
				}
				;
				parallel_ir_array_temp__9_62_1 = labels.ARRAYELEM(_dist_parfor_95_index2,((_dist_parfor_94_index1) - (__hpat_loop_start_100)) + (1));
				SSAValue15 = -(parallel_ir_array_temp__9_62_1);
				parallel_ir_array_temp_SSAValue53_67_1 = _dist_gemm_tmp3_99;
				SSAValue16 = (SSAValue15) * (parallel_ir_array_temp_SSAValue53_67_1);
				SSAValue17 = exp(SSAValue16);
				SSAValue19 = (1.0) + (SSAValue17);
				SSAValue20 = (1.0) / (SSAValue19);
				SSAValue21 = (SSAValue20) - (1.0);
				parallel_ir_array_temp__9_88_1 = labels.ARRAYELEM(_dist_parfor_95_index2,((_dist_parfor_94_index1) - (__hpat_loop_start_100)) + (1));
				SSAValue22 = (SSAValue21) * (parallel_ir_array_temp__9_88_1);
				_dist_gemm_tmp2_105 = SSAValue22;
				for(int64_t _dist_loop_103_index = 1; _dist_loop_103_index <= __hpat_h5_dim_size_1_2; _dist_loop_103_index += 1)
				{
					;
					_dist_gemm_tmp1_104 = points.ARRAYELEM(_dist_loop_103_index,((_dist_parfor_94_index1) - (__hpat_loop_start_100)) + (1));
					_dist_gemm_tmp3_106 = SSAValue13.ARRAYELEM(_dist_parfor_95_index2,_dist_loop_103_index);
					_dist_gemm_tmp3_106 = (_dist_gemm_tmp3_106) + ((_dist_gemm_tmp1_104) * (_dist_gemm_tmp2_105));
					SSAValue13.ARRAYELEM(_dist_parfor_95_index2,_dist_loop_103_index) = _dist_gemm_tmp3_106;
				}
				;
			}
		}
		;
		__hpat_reduce_111 = j2c_array<double>::new_j2c_array_2d(NULL, 1, __hpat_h5_dim_size_1_2);
		MPI_Allreduce(SSAValue13.data, __hpat_reduce_111.data, (1) * (__hpat_h5_dim_size_1_2), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
		SSAValue13 = __hpat_reduce_111;
		parallel_ir_new_array_name_91_1 = j2c_array<double>::new_j2c_array_2d(NULL, 1, __hpat_h5_dim_size_1_2);
		for ( parfor_index_2_91 = 1; parfor_index_2_91 <= (int64_t)__hpat_h5_dim_size_1_2; parfor_index_2_91 += 1)
		{
			for ( parfor_index_1_91 = 1; parfor_index_1_91 <= (int64_t)1; parfor_index_1_91 += 1)
			{
				;
				parallel_ir_array_temp__10_92_1 = w.ARRAYELEM(parfor_index_1_91,parfor_index_2_91);
				parallel_ir_array_temp_SSAValue61_93_1 = SSAValue13.ARRAYELEM(parfor_index_1_91,parfor_index_2_91);
				SSAValue1 = (parallel_ir_array_temp__10_92_1) - (parallel_ir_array_temp_SSAValue61_93_1);
				parallel_ir_new_array_name_91_1.ARRAYELEM(parfor_index_1_91,parfor_index_2_91) = SSAValue1;
			}
		}
		;
		w = parallel_ir_new_array_name_91_1;
	}
	;
	0; if(__hpat_node_id==0) printf("exec time %lf\n", MPI_Wtime()-__hpat_t1);;
	*ret0 = w;
	return;

}


extern "C" void _pplogistic_regressionp271_(int run_where, int64_t iterations, ASCIIString&  file_name ,  j2c_array< double > ** __restrict ret0 , bool genMain = true)
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
		mainFileData << iterations << std::endl;
		mainFileData << file_name << std::endl;
		mainFileData.close();
		std::ofstream mainFile(newMain.str());
		mainFile << "#include \"" << __FILE__ << "\"" << std::endl;
		mainFile << "int main(int argc, char *argv[]) {" << std::endl;
		mainFile << "    MPI_Init(&argc, &argv);" << std::endl;
		mainFile << "    std::ifstream mainFileData(\"" << newMainData.str() << "\", std::ios::in | std::ios::binary);" << std::endl;
		mainFile << "    int runwhere;" << std::endl;
		mainFile << "    mainFileData >> runwhere;" << std::endl;
		mainFile << "     j2c_array< double > * ret0;" << std::endl;
		mainFile << "    int64_t iterations;" << std::endl;
		mainFile << "    mainFileData >> iterations;" << std::endl;
		mainFile << "    ASCIIString file_name;" << std::endl;
		mainFile << "    mainFileData >> file_name;" << std::endl;
		mainFile << "    _pplogistic_regressionp271_(runwhere, iterations, file_name, &ret0, false);" << std::endl;
		mainFile << "    MPI_Finalize();" << std::endl;
		mainFile << "    return 0;" << std::endl;
		mainFile << "}" << std::endl;
		mainFile.close();
		std::ofstream mainFileSh(newMainSh.str());
		mainFileSh << "#!/bin/sh" << std::endl;
		mainFileSh << "mpiicpc -O3 -std=c++11 -I/usr/local/hdf5/include  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << " -L/usr/local/hdf5/lib -lhdf5  -lm " << std::endl;
		mainFileSh.close();
	}

	*ret0 = new  j2c_array< double > ();

	pplogistic_regressionp271(iterations, file_name, *ret0);
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
		case 2:
			a = new  j2c_array< double > ((double*)data, ndim, dims);
			break;
		default:
			fprintf(stderr, "j2c_array_new called with invalid key %d", key);
			assert(false);
			break;
	}
	return a;
}
