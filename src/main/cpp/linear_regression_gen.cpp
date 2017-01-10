#include <random>
#include <mpi.h>
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
} ppplinear_regressionp271;
void pplinear_regressionp271(int64_t iterations, int64_t N,  j2c_array< double >  * __restrict ret0)
{
	ppplinear_regressionp271 pselfp;
	j2c_array< double >  w;
	double alphaN;
	int64_t ptempp;
	j2c_array< double >  _pa_rand_gen_arr;
	j2c_array< double >  pp_pa_rand_gen_arr_14p276;
	double SSAValue40pp2;
	j2c_array< double >  SSAValue0;
	j2c_array< double >  SSAValue1;
	int64_t SSAValue2;
	j2c_array< double >  SSAValue4;
	double SSAValue5;
	bool SSAValue6;
	int64_t SSAValue7;
	int64_t SSAValue8;
	bool SSAValue9;
	bool SSAValue10;
	j2c_array< double >  SSAValue11;
	j2c_array< double >  SSAValue12;
	int64_t parfor_index_1_29;
	int64_t parfor_index_2_29;
	double SSAValue13;
	int64_t parfor_index_1_33;
	int64_t parfor_index_2_33;
	double SSAValue14;
	int64_t parfor_index_1_37;
	int64_t parfor_index_2_37;
	double parallel_ir_array_temp_SSAValue55_42_1;
	double parallel_ir_array_temp__12_43_1;
	int64_t parfor_index_1_41;
	int64_t parfor_index_2_41;
	double SSAValue15;
	j2c_array< double >  parallel_ir_new_array_name_41_1;
	double SSAValue16;
	double parallel_ir_array_temp__8_51_1;
	double parallel_ir_array_temp_SSAValue59_52_1;
	int64_t parfor_index_1_50;
	int64_t parfor_index_2_50;
	double SSAValue3;
	j2c_array< double >  parallel_ir_new_array_name_50_1;
	int64_t _dist_parfor_53_index1;
	int64_t _dist_parfor_54_index2;
	int64_t _dist_loop_55_index;
	double _dist_gemm_tmp1_56;
	double _dist_gemm_tmp2_57;
	double _dist_gemm_tmp3_58;
	int64_t _dist_parfor_60_index1;
	int64_t _dist_parfor_61_index2;
	int64_t _dist_loop_62_index;
	double _dist_gemm_tmp1_63;
	double _dist_gemm_tmp2_64;
	double _dist_gemm_tmp3_65;
	int32_t __hpat_num_pes;
	int32_t __hpat_node_id;
	int64_t __hpat_loop_start_29;
	int64_t __hpat_loop_end_29;
	int64_t __hpat_loop_div_29;
	int64_t __hpat_loop_start_33;
	int64_t __hpat_loop_end_33;
	int64_t __hpat_loop_div_33;
	int64_t __hpat_loop_start_59;
	int64_t __hpat_loop_end_59;
	int64_t __hpat_loop_div_59;
	j2c_array< double >  __hpat_reduce_67;
	std::random_device cgen_rand_device;
	std::uniform_real_distribution<double> cgen_distribution(0.0,1.0);
	std::normal_distribution<double> cgen_n_distribution(0.0,1.0);
	//std::default_random_engine cgen_rand_generator(cgen_rand_device());
	std::default_random_engine cgen_rand_generator(0);
	;;
	MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
	MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
	SSAValue40pp2 = (double)0;
	SSAValue5 = (double)N;
	alphaN = (0.01) / (SSAValue5);
	SSAValue6 = (1) <= (iterations);
	SSAValue7 = (1) - (1);
	SSAValue2 = (SSAValue6) ? (iterations) : (SSAValue7);
	SSAValue8 = (SSAValue2) + (1);
	ptempp = 1;
	__hpat_loop_div_29 = (N) / (__hpat_num_pes);
	__hpat_loop_start_29 = ((__hpat_node_id) * (__hpat_loop_div_29)) + (1);
	__hpat_loop_end_29 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_29);
	_pa_rand_gen_arr = j2c_array<double>::new_j2c_array_2d(NULL, 4, N);
	#pragma simd
	for ( parfor_index_2_29 = __hpat_loop_start_29; parfor_index_2_29 <= (int64_t)__hpat_loop_end_29; parfor_index_2_29 += 1)
	{
		for ( parfor_index_1_29 = 1; parfor_index_1_29 <= (int64_t)4; parfor_index_1_29 += 1)
		{
			;
			SSAValue13 = cgen_distribution(cgen_rand_generator);
			;
			_pa_rand_gen_arr.ARRAYELEM(parfor_index_1_29,((parfor_index_2_29) - (__hpat_loop_start_29)) + (1)) = SSAValue13;
		}
	}
	;
	__hpat_loop_div_33 = (N) / (__hpat_num_pes);
	__hpat_loop_start_33 = ((__hpat_node_id) * (__hpat_loop_div_33)) + (1);
	__hpat_loop_end_33 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_33);
	pp_pa_rand_gen_arr_14p276 = j2c_array<double>::new_j2c_array_2d(NULL, 10, N);
	#pragma simd
	for ( parfor_index_2_33 = __hpat_loop_start_33; parfor_index_2_33 <= (int64_t)__hpat_loop_end_33; parfor_index_2_33 += 1)
	{
		for ( parfor_index_1_33 = 1; parfor_index_1_33 <= (int64_t)10; parfor_index_1_33 += 1)
		{
			;
			SSAValue14 = cgen_distribution(cgen_rand_generator);
			;
			pp_pa_rand_gen_arr_14p276.ARRAYELEM(parfor_index_1_33,((parfor_index_2_33) - (__hpat_loop_start_33)) + (1)) = SSAValue14;
		}
	}
	;
	SSAValue4 = j2c_array<double>::new_j2c_array_2d(NULL, 4, 10);
	for ( parfor_index_2_37 = 1; parfor_index_2_37 <= (int64_t)10; parfor_index_2_37 += 1)
	{
		for ( parfor_index_1_37 = 1; parfor_index_1_37 <= (int64_t)4; parfor_index_1_37 += 1)
		{
			;
			SSAValue4.ARRAYELEM(parfor_index_1_37,parfor_index_2_37) = SSAValue40pp2;
		}
	}
	;
	w = SSAValue4;
	label37 : ;
	SSAValue9 = (ptempp == SSAValue8);
	SSAValue10 = !(SSAValue9);
	if (!(SSAValue10)) goto label71;
	ptempp = (ptempp) + (1);
	__hpat_loop_div_59 = (N) / (__hpat_num_pes);
	__hpat_loop_start_59 = ((__hpat_node_id) * (__hpat_loop_div_59)) + (1);
	__hpat_loop_end_59 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_59);
	SSAValue12 = j2c_array<double>::new_j2c_array_2d(NULL, 4, N);
	memset(SSAValue12.data, 0, sizeof(double)*(10) * (4));
	#pragma simd
	for ( _dist_parfor_53_index1 = __hpat_loop_start_59; _dist_parfor_53_index1 <= (int64_t)__hpat_loop_end_59; _dist_parfor_53_index1 += 1)
	{
		for ( _dist_parfor_54_index2 = 1; _dist_parfor_54_index2 <= (int64_t)4; _dist_parfor_54_index2 += 1)
		{
			;
			_dist_gemm_tmp3_58 = 0;
			for(int64_t _dist_loop_55_index = 1; _dist_loop_55_index <= 10; _dist_loop_55_index += 1)
			{
				;
				_dist_gemm_tmp1_56 = w.ARRAYELEM(_dist_parfor_54_index2,_dist_loop_55_index);
				_dist_gemm_tmp2_57 = pp_pa_rand_gen_arr_14p276.ARRAYELEM(_dist_loop_55_index,((_dist_parfor_53_index1) - (__hpat_loop_start_59)) + (1));
				_dist_gemm_tmp3_58 = (_dist_gemm_tmp3_58) + ((_dist_gemm_tmp1_56) * (_dist_gemm_tmp2_57));
			}
			;
			parallel_ir_array_temp_SSAValue55_42_1 = _dist_gemm_tmp3_58;
			parallel_ir_array_temp__12_43_1 = _pa_rand_gen_arr.ARRAYELEM(_dist_parfor_54_index2,((_dist_parfor_53_index1) - (__hpat_loop_start_59)) + (1));
			SSAValue15 = (parallel_ir_array_temp_SSAValue55_42_1) - (parallel_ir_array_temp__12_43_1);
			SSAValue16 = (alphaN) * (SSAValue15);
			_dist_gemm_tmp2_64 = SSAValue16;
			for(int64_t _dist_loop_62_index = 1; _dist_loop_62_index <= 10; _dist_loop_62_index += 1)
			{
				;
				_dist_gemm_tmp1_63 = pp_pa_rand_gen_arr_14p276.ARRAYELEM(_dist_loop_62_index,((_dist_parfor_53_index1) - (__hpat_loop_start_59)) + (1));
				_dist_gemm_tmp3_65 = SSAValue12.ARRAYELEM(_dist_parfor_54_index2,_dist_loop_62_index);
				_dist_gemm_tmp3_65 = (_dist_gemm_tmp3_65) + ((_dist_gemm_tmp1_63) * (_dist_gemm_tmp2_64));
				SSAValue12.ARRAYELEM(_dist_parfor_54_index2,_dist_loop_62_index) = _dist_gemm_tmp3_65;
			}
			;
		}
	}
	;
	__hpat_reduce_67 = j2c_array<double>::new_j2c_array_2d(NULL, 4, 10);
	MPI_Allreduce(SSAValue12.data, __hpat_reduce_67.data, (4) * (10), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
	SSAValue12 = __hpat_reduce_67;
	parallel_ir_new_array_name_50_1 = j2c_array<double>::new_j2c_array_2d(NULL, 4, 10);
	for ( parfor_index_2_50 = 1; parfor_index_2_50 <= (int64_t)10; parfor_index_2_50 += 1)
	{
		for ( parfor_index_1_50 = 1; parfor_index_1_50 <= (int64_t)4; parfor_index_1_50 += 1)
		{
			;
			parallel_ir_array_temp__8_51_1 = w.ARRAYELEM(parfor_index_1_50,parfor_index_2_50);
			parallel_ir_array_temp_SSAValue59_52_1 = SSAValue12.ARRAYELEM(parfor_index_1_50,parfor_index_2_50);
			SSAValue3 = (parallel_ir_array_temp__8_51_1) - (parallel_ir_array_temp_SSAValue59_52_1);
			parallel_ir_new_array_name_50_1.ARRAYELEM(parfor_index_1_50,parfor_index_2_50) = SSAValue3;
		}
	}
	;
	w = parallel_ir_new_array_name_50_1;
	goto label37;
	label71 : ;
	*ret0 = w;
	return;

}


extern "C" void _pplinear_regressionp271_(int run_where, int64_t iterations, int64_t N ,  j2c_array< double > ** __restrict ret0 , bool genMain = true)
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
		mainFileData << N << std::endl;
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
		mainFile << "    int64_t N;" << std::endl;
		mainFile << "    mainFileData >> N;" << std::endl;
		mainFile << "    _pplinear_regressionp271_(runwhere, iterations, N, &ret0, false);" << std::endl;
		mainFile << "    MPI_Finalize();" << std::endl;
		mainFile << "    return 0;" << std::endl;
		mainFile << "}" << std::endl;
		mainFile.close();
		std::ofstream mainFileSh(newMainSh.str());
		mainFileSh << "#!/bin/sh" << std::endl;
		mainFileSh << "mpiicpc -O3 -std=c++11 -I/usr/local/hdf5/include  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << "  -lm " << std::endl;
		mainFileSh.close();
	}

	*ret0 = new  j2c_array< double > ();

	pplinear_regressionp271(iterations, N, *ret0);
}


extern "C"
void *j2c_array_new(int key, void*data, unsigned ndim, int64_t *dims)
{
	void *a = NULL;
	switch(key)
	{
		case 1:
			a = new  j2c_array< double > ((double*)data, ndim, dims);
			break;
		default:
			fprintf(stderr, "j2c_array_new called with invalid key %d", key);
			assert(false);
			break;
	}
	return a;
}
