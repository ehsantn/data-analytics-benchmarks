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
} ppplogistic_regressionp271;
typedef struct
{
	int64_t f0;
	int64_t f1;
} TupleInt64Int64;
void pplogistic_regressionp271(int64_t iterations, int64_t N,  j2c_array< double >  * __restrict ret0)
{
	ppplogistic_regressionp271 pselfp;
	j2c_array< double >  w;
	int64_t ptempp;
	j2c_array< double >  _pa_rand_gen_arr;
	j2c_array< double >  pp_pa_rand_gen_arr_12p280;
	j2c_array< double >  pp_pa_rand_gen_arr_14p282;
	j2c_array< double >  SSAValue17;
	j2c_array< double >  SSAValue23;
	int64_t SSAValue0;
	TupleInt64Int64 SSAValue2;
	bool SSAValue3;
	int64_t SSAValue4;
	int64_t SSAValue5;
	bool SSAValue6;
	bool SSAValue7;
	j2c_array< double >  SSAValue8;
	j2c_array< double >  SSAValue9;
	int64_t parfor_index_1_59;
	int64_t parfor_index_2_59;
	double SSAValue10;
	int64_t parfor_index_1_63;
	int64_t parfor_index_2_63;
	double SSAValue11;
	int64_t parfor_index_1_67;
	double SSAValue13;
	double SSAValue14;
	double SSAValue15;
	double parallel_ir_array_temp__10_80_1;
	int64_t parfor_index_1_79;
	int64_t parfor_index_2_79;
	double SSAValue16;
	j2c_array< double >  parallel_ir_new_array_name_79_1;
	double parallel_ir_array_temp_SSAValue47_85_1;
	double SSAValue18;
	double SSAValue19;
	double SSAValue20;
	double SSAValue21;
	double SSAValue22;
	double parallel_ir_array_temp__10_106_1;
	double SSAValue1;
	double parallel_ir_array_temp__7_110_1;
	double parallel_ir_array_temp_SSAValue55_111_1;
	int64_t parfor_index_1_109;
	int64_t parfor_index_2_109;
	double SSAValue12;
	j2c_array< double >  parallel_ir_new_array_name_109_1;
	int64_t _dist_parfor_112_index1;
	int64_t _dist_parfor_113_index2;
	int64_t _dist_loop_114_index;
	double _dist_gemm_tmp1_115;
	double _dist_gemm_tmp2_116;
	double _dist_gemm_tmp3_117;
	int64_t _dist_parfor_119_index1;
	int64_t _dist_parfor_120_index2;
	int64_t _dist_loop_121_index;
	double _dist_gemm_tmp1_122;
	double _dist_gemm_tmp2_123;
	double _dist_gemm_tmp3_124;
	int32_t __hpat_num_pes;
	int32_t __hpat_node_id;
	int64_t __hpat_loop_start_59;
	int64_t __hpat_loop_end_59;
	int64_t __hpat_loop_div_59;
	int64_t __hpat_loop_start_63;
	int64_t __hpat_loop_end_63;
	int64_t __hpat_loop_div_63;
	int64_t __hpat_bcast_size_75;
	int64_t __hpat_loop_start_118;
	int64_t __hpat_loop_end_118;
	int64_t __hpat_loop_div_118;
	j2c_array< double >  __hpat_reduce_126;
	std::random_device cgen_rand_device;
	std::uniform_real_distribution<double> cgen_distribution(0.0,1.0);
	std::normal_distribution<double> cgen_n_distribution(0.0,1.0);
	std::default_random_engine cgen_rand_generator(cgen_rand_device());
	;;
	MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
	MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
	SSAValue2 = TupleInt64Int64{1, 10};
	SSAValue3 = (1) <= (iterations);
	SSAValue4 = (1) - (1);
	SSAValue0 = (SSAValue3) ? (iterations) : (SSAValue4);
	SSAValue5 = (SSAValue0) + (1);
	__hpat_loop_div_59 = (N) / (__hpat_num_pes);
	__hpat_loop_start_59 = ((__hpat_node_id) * (__hpat_loop_div_59)) + (1);
	__hpat_loop_end_59 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_59);
	_pa_rand_gen_arr = j2c_array<double>::new_j2c_array_2d(NULL, 1, N);
	#pragma simd
	for ( parfor_index_2_59 = __hpat_loop_start_59; parfor_index_2_59 <= (int64_t)__hpat_loop_end_59; parfor_index_2_59 += 1)
	{
		for ( parfor_index_1_59 = 1; parfor_index_1_59 <= (int64_t)1; parfor_index_1_59 += 1)
		{
			;
			SSAValue10 = cgen_distribution(cgen_rand_generator);
			;
			_pa_rand_gen_arr.ARRAYELEM(parfor_index_1_59,((parfor_index_2_59) - (__hpat_loop_start_59)) + (1)) = SSAValue10;
		}
	}
	;
	__hpat_loop_div_63 = (N) / (__hpat_num_pes);
	__hpat_loop_start_63 = ((__hpat_node_id) * (__hpat_loop_div_63)) + (1);
	__hpat_loop_end_63 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_63);
	pp_pa_rand_gen_arr_12p280 = j2c_array<double>::new_j2c_array_2d(NULL, 10, N);
	#pragma simd
	for ( parfor_index_2_63 = __hpat_loop_start_63; parfor_index_2_63 <= (int64_t)__hpat_loop_end_63; parfor_index_2_63 += 1)
	{
		for ( parfor_index_1_63 = 1; parfor_index_1_63 <= (int64_t)10; parfor_index_1_63 += 1)
		{
			;
			SSAValue11 = cgen_distribution(cgen_rand_generator);
			;
			pp_pa_rand_gen_arr_12p280.ARRAYELEM(parfor_index_1_63,((parfor_index_2_63) - (__hpat_loop_start_63)) + (1)) = SSAValue11;
		}
	}
	;
	pp_pa_rand_gen_arr_14p282 = j2c_array<double>::new_j2c_array_1d(NULL, 10);
	if ((__hpat_node_id == 0))
	{
		;
		for ( parfor_index_1_67 = 1; parfor_index_1_67 <= (int64_t)10; parfor_index_1_67 += 1)
		{
			;
			SSAValue13 = cgen_distribution(cgen_rand_generator);
			;
			SSAValue14 = (2.0) * (SSAValue13);
			SSAValue15 = (SSAValue14) - (1.0);
			pp_pa_rand_gen_arr_14p282.ARRAYELEM(parfor_index_1_67) = SSAValue15;
		}
		;
	}
	;
	__hpat_bcast_size_75 = 10;
	MPI_Bcast(pp_pa_rand_gen_arr_14p282.data, __hpat_bcast_size_75, MPI_DOUBLE, 0, MPI_COMM_WORLD);;
	w = pp_pa_rand_gen_arr_14p282.reshape(SSAValue2.f0,SSAValue2.f1);
	;
	0; double __hpat_t1 = MPI_Wtime();
	ptempp = 1;
	while (1)
	{
		;
		SSAValue6 = (ptempp == SSAValue5);
		SSAValue7 = !(SSAValue6);
		if (!(SSAValue7)) break;
		ptempp = (ptempp) + (1);
		__hpat_loop_div_118 = (N) / (__hpat_num_pes);
		__hpat_loop_start_118 = ((__hpat_node_id) * (__hpat_loop_div_118)) + (1);
		__hpat_loop_end_118 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_118);
		SSAValue9 = j2c_array<double>::new_j2c_array_2d(NULL, 1, N);
		memset(SSAValue9.data, 0, sizeof(double)*(10) * (1));
		#pragma simd
		for ( _dist_parfor_112_index1 = __hpat_loop_start_118; _dist_parfor_112_index1 <= (int64_t)__hpat_loop_end_118; _dist_parfor_112_index1 += 1)
		{
			for ( _dist_parfor_113_index2 = 1; _dist_parfor_113_index2 <= (int64_t)1; _dist_parfor_113_index2 += 1)
			{
				;
				_dist_gemm_tmp3_117 = 0;
				for(int64_t _dist_loop_114_index = 1; _dist_loop_114_index <= 10; _dist_loop_114_index += 1)
				{
					;
					_dist_gemm_tmp1_115 = w.ARRAYELEM(_dist_parfor_113_index2,_dist_loop_114_index);
					_dist_gemm_tmp2_116 = pp_pa_rand_gen_arr_12p280.ARRAYELEM(_dist_loop_114_index,((_dist_parfor_112_index1) - (__hpat_loop_start_118)) + (1));
					_dist_gemm_tmp3_117 = (_dist_gemm_tmp3_117) + ((_dist_gemm_tmp1_115) * (_dist_gemm_tmp2_116));
				}
				;
				parallel_ir_array_temp__10_80_1 = _pa_rand_gen_arr.ARRAYELEM(_dist_parfor_113_index2,((_dist_parfor_112_index1) - (__hpat_loop_start_118)) + (1));
				SSAValue16 = -(parallel_ir_array_temp__10_80_1);
				parallel_ir_array_temp_SSAValue47_85_1 = _dist_gemm_tmp3_117;
				SSAValue18 = (SSAValue16) * (parallel_ir_array_temp_SSAValue47_85_1);
				SSAValue19 = exp(SSAValue18);
				SSAValue20 = (1.0) + (SSAValue19);
				SSAValue21 = (1.0) / (SSAValue20);
				SSAValue22 = (SSAValue21) - (1.0);
				parallel_ir_array_temp__10_106_1 = _pa_rand_gen_arr.ARRAYELEM(_dist_parfor_113_index2,((_dist_parfor_112_index1) - (__hpat_loop_start_118)) + (1));
				SSAValue1 = (SSAValue22) * (parallel_ir_array_temp__10_106_1);
				_dist_gemm_tmp2_123 = SSAValue1;
				for(int64_t _dist_loop_121_index = 1; _dist_loop_121_index <= 10; _dist_loop_121_index += 1)
				{
					;
					_dist_gemm_tmp1_122 = pp_pa_rand_gen_arr_12p280.ARRAYELEM(_dist_loop_121_index,((_dist_parfor_112_index1) - (__hpat_loop_start_118)) + (1));
					_dist_gemm_tmp3_124 = SSAValue9.ARRAYELEM(_dist_parfor_113_index2,_dist_loop_121_index);
					_dist_gemm_tmp3_124 = (_dist_gemm_tmp3_124) + ((_dist_gemm_tmp1_122) * (_dist_gemm_tmp2_123));
					SSAValue9.ARRAYELEM(_dist_parfor_113_index2,_dist_loop_121_index) = _dist_gemm_tmp3_124;
				}
				;
			}
		}
		;
		__hpat_reduce_126 = j2c_array<double>::new_j2c_array_2d(NULL, 1, 10);
		MPI_Allreduce(SSAValue9.data, __hpat_reduce_126.data, (1) * (10), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
		SSAValue9 = __hpat_reduce_126;
		parallel_ir_new_array_name_109_1 = j2c_array<double>::new_j2c_array_2d(NULL, 1, 10);
		for ( parfor_index_2_109 = 1; parfor_index_2_109 <= (int64_t)10; parfor_index_2_109 += 1)
		{
			for ( parfor_index_1_109 = 1; parfor_index_1_109 <= (int64_t)1; parfor_index_1_109 += 1)
			{
				;
				parallel_ir_array_temp__7_110_1 = w.ARRAYELEM(parfor_index_1_109,parfor_index_2_109);
				parallel_ir_array_temp_SSAValue55_111_1 = SSAValue9.ARRAYELEM(parfor_index_1_109,parfor_index_2_109);
				SSAValue12 = (parallel_ir_array_temp__7_110_1) - (parallel_ir_array_temp_SSAValue55_111_1);
				parallel_ir_new_array_name_109_1.ARRAYELEM(parfor_index_1_109,parfor_index_2_109) = SSAValue12;
			}
		}
		;
		w = parallel_ir_new_array_name_109_1;
	}
	;
	0; if(__hpat_node_id==0) printf("exec time %lf\n", MPI_Wtime()-__hpat_t1);;
	*ret0 = w;
	return;

}


extern "C" void _pplogistic_regressionp271_(int run_where, int64_t iterations, int64_t N ,  j2c_array< double > ** __restrict ret0 , bool genMain = true)
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
		mainFile << "    _pplogistic_regressionp271_(runwhere, iterations, N, &ret0, false);" << std::endl;
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

	pplogistic_regressionp271(iterations, N, *ret0);
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
