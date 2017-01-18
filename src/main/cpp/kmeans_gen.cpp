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
} pppkmeansp271;
typedef struct
{
	double val;
	int64_t idx;
} ParallelAcceleratorAPILibMTFloat64;
void ppkmeansp271(int64_t numCenter, int64_t iterNum, int64_t N,  j2c_array< double >  * __restrict ret0)
{
	pppkmeansp271 pselfp;
	j2c_array< double >  centroids;
	int64_t ptempp;
	j2c_array< double >  _pa_rand_gen_arr;
	j2c_array< double >  pp_pa_rand_gen_arr_17p285;
	int64_t SSAValue29;
	bool SSAValue14;
	int64_t SSAValue15;
	int64_t SSAValue27;
	bool SSAValue0;
	bool SSAValue1;
	j2c_array< int64_t >  SSAValue3;
	j2c_array< double >  SSAValue4;
	int64_t parfor_index_1_17;
	int64_t parfor_index_2_17;
	double SSAValue5;
	int64_t parfor_index_1_21;
	int64_t parfor_index_2_21;
	double SSAValue6;
	int64_t parfor_index_1_25;
	j2c_array< double >  SSAValue7;
	int64_t parfor_index_1_30;
	double SSAValue10;
	double parallel_ir_temp__8_parallel_ir_range__37Int64_skip__39Int64_1_1_1;
	double parallel_ir_temp__6_parallel_ir_range__43Int64_skip__45Int64_1_1_1;
	int64_t parfor_index_1_51;
	double SSAValue11;
	int32_t SSAValue12;
	double SSAValue13;
	double parallel_ir_array_temp_SSAValue3_68_1;
	double parallel_ir_reduction_output_66;
	ParallelAcceleratorAPILibMTFloat64 m;
	int64_t SSAValue16;
	int64_t ppip300;
	ParallelAcceleratorAPILibMTFloat64 parallel_ir_reduction_input_77_1;
	ParallelAcceleratorAPILibMTFloat64 pptemp_neutral_valp301;
	double SSAValue17;
	int64_t SSAValue18;
	double SSAValue19;
	double SSAValue20;
	bool SSAValue21;
	double SSAValue22;
	double SSAValue23;
	double SSAValue24;
	bool SSAValue25;
	int64_t parfor_index_1_82;
	int64_t parfor_index_2_82;
	double SSAValue28;
	double SSAValue2;
	int64_t parallel_ir_array_temp__7_92_1;
	int64_t parfor_index_1_91;
	bool SSAValue8;
	bool parallel_ir_array_temp__21_94_1;
	double parallel_ir_array_temp__9_97_1;
	double parallel_ir_reduction_output_95;
	int64_t parallel_ir_array_temp__7_100_1;
	bool SSAValue9;
	int64_t SSAValue26;
	int64_t parallel_ir_array_temp_SSAValue6_109_1;
	int64_t parallel_ir_reduction_output_107;
	j2c_array< double >  HPAT__99_82;
	j2c_array< int64_t >  HPAT__105_82;
	j2c_array< double >  HPAT_SSAValue28_82;
	j2c_array< double >  HPAT_SSAValue2_82;
	int32_t __hpat_num_pes;
	int32_t __hpat_node_id;
	int64_t __hpat_dist_arr_start_94;
	int64_t __hpat_dist_arr_div_94;
	int64_t __hpat_dist_arr_count_94;
	int64_t __hpat_loop_start_17;
	int64_t __hpat_loop_end_17;
	int64_t __hpat_loop_div_17;
	int64_t __hpat_bcast_size_107;
	int64_t __hpat_loop_start_25;
	int64_t __hpat_loop_end_25;
	int64_t __hpat_loop_div_25;
	j2c_array< double >  __hpat_reduce_95;
	j2c_array< int64_t >  __hpat_reduce_96;
	std::random_device cgen_rand_device;
	std::uniform_real_distribution<double> cgen_distribution(0.0,1.0);
	std::normal_distribution<double> cgen_n_distribution(0.0,1.0);
	std::default_random_engine cgen_rand_generator(cgen_rand_device());
	;;
	MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
	MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
	SSAValue14 = (1) <= (iterNum);
	SSAValue15 = (1) - (1);
	SSAValue29 = (SSAValue14) ? (iterNum) : (SSAValue15);
	SSAValue27 = (SSAValue29) + (1);
	__hpat_loop_div_17 = (N) / (__hpat_num_pes);
	__hpat_loop_start_17 = ((__hpat_node_id) * (__hpat_loop_div_17)) + (1);
	__hpat_loop_end_17 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_17);
	__hpat_dist_arr_div_94 = (N) / (__hpat_num_pes);
	__hpat_dist_arr_start_94 = (__hpat_node_id) * (__hpat_dist_arr_div_94);
	__hpat_dist_arr_count_94 = ((__hpat_node_id==__hpat_num_pes-1) ? N-__hpat_node_id*__hpat_dist_arr_div_94 : __hpat_dist_arr_div_94);
	_pa_rand_gen_arr = j2c_array<double>::new_j2c_array_2d(NULL, 10, __hpat_dist_arr_count_94);
	#pragma simd
	for ( parfor_index_2_17 = __hpat_loop_start_17; parfor_index_2_17 <= (int64_t)__hpat_loop_end_17; parfor_index_2_17 += 1)
	{
		for ( parfor_index_1_17 = 1; parfor_index_1_17 <= (int64_t)10; parfor_index_1_17 += 1)
		{
			;
			SSAValue5 = cgen_distribution(cgen_rand_generator);
			;
			_pa_rand_gen_arr.ARRAYELEM(parfor_index_1_17,((parfor_index_2_17) - (__hpat_loop_start_17)) + (1)) = SSAValue5;
		}
	}
	;
	pp_pa_rand_gen_arr_17p285 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
	if ((__hpat_node_id == 0))
	{
		;
		for ( parfor_index_2_21 = 1; parfor_index_2_21 <= (int64_t)numCenter; parfor_index_2_21 += 1)
		{
			for ( parfor_index_1_21 = 1; parfor_index_1_21 <= (int64_t)10; parfor_index_1_21 += 1)
			{
				;
				SSAValue6 = cgen_distribution(cgen_rand_generator);
				;
				pp_pa_rand_gen_arr_17p285.ARRAYELEM(parfor_index_1_21,parfor_index_2_21) = SSAValue6;
			}
		}
		;
	}
	;
	__hpat_bcast_size_107 = (10) * (numCenter);
	MPI_Bcast(pp_pa_rand_gen_arr_17p285.data, __hpat_bcast_size_107, MPI_DOUBLE, 0, MPI_COMM_WORLD);;
	centroids = pp_pa_rand_gen_arr_17p285;
	0; double __hpat_t1 = MPI_Wtime();
	ptempp = 1;
	while (1)
	{
		;
		SSAValue0 = (ptempp == SSAValue27);
		SSAValue1 = !(SSAValue0);
		if (!(SSAValue1)) break;
		ptempp = (ptempp) + (1);
		SSAValue4 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
		__hpat_loop_div_25 = (N) / (__hpat_num_pes);
		__hpat_loop_start_25 = ((__hpat_node_id) * (__hpat_loop_div_25)) + (1);
		__hpat_loop_end_25 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_25);
		HPAT__99_82 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
		HPAT__105_82 = j2c_array<int64_t>::new_j2c_array_2d(NULL, 10, numCenter);
		HPAT_SSAValue28_82 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
		HPAT_SSAValue2_82 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
		for(int64_t parfor_index_2_82 = 1; parfor_index_2_82 <= numCenter; parfor_index_2_82 += 1)
		{
			;
			for(int64_t parfor_index_1_82 = 1; parfor_index_1_82 <= 10; parfor_index_1_82 += 1)
			{
				;
				HPAT__99_82.ARRAYELEM(parfor_index_1_82, parfor_index_2_82) = 0.0;
				HPAT__105_82.ARRAYELEM(parfor_index_1_82, parfor_index_2_82) = 0;
			}
			;
		}
		;
		#pragma simd
		for ( parfor_index_1_25 = __hpat_loop_start_25; parfor_index_1_25 <= (int64_t)__hpat_loop_end_25; parfor_index_1_25 += 1)
		{
			;
			m = ParallelAcceleratorAPILibMTFloat64{DBL_MAX, 1};
			for ( parfor_index_1_30 = 1; parfor_index_1_30 <= (int64_t)numCenter; parfor_index_1_30 += 1)
			{
				;
				parallel_ir_reduction_output_66 = 0.0;
				for ( parfor_index_1_51 = 1; parfor_index_1_51 <= (int64_t)10; parfor_index_1_51 += 1)
				{
					;
					parallel_ir_temp__8_parallel_ir_range__37Int64_skip__39Int64_1_1_1 = _pa_rand_gen_arr.ARRAYELEM(parfor_index_1_51,((parfor_index_1_25) - (__hpat_loop_start_25)) + (1));
					parallel_ir_temp__6_parallel_ir_range__43Int64_skip__45Int64_1_1_1 = centroids.ARRAYELEM(parfor_index_1_51,parfor_index_1_30);
					SSAValue11 = (parallel_ir_temp__8_parallel_ir_range__37Int64_skip__39Int64_1_1_1) - (parallel_ir_temp__6_parallel_ir_range__43Int64_skip__45Int64_1_1_1);
					SSAValue12 = (int32_t)(2);
					SSAValue13 = pow(SSAValue11, SSAValue12);
					parallel_ir_reduction_output_66 = (parallel_ir_reduction_output_66+SSAValue13);
				}
				;
				SSAValue10 = sqrt(parallel_ir_reduction_output_66);;
				SSAValue23 = SSAValue10;
				SSAValue24 = m.val;
				SSAValue25 = (SSAValue23) < (SSAValue24);
				if (SSAValue25)
				{
					;
					SSAValue22 = SSAValue10;
					m.val = SSAValue22;
					m.idx = parfor_index_1_30;
				}
				;
			}
			;
			SSAValue16 = m.idx;
			for ( parfor_index_2_82 = 1; parfor_index_2_82 <= (int64_t)numCenter; parfor_index_2_82 += 1)
			{
				for ( parfor_index_1_82 = 1; parfor_index_1_82 <= (int64_t)10; parfor_index_1_82 += 1)
				{
					;
					parallel_ir_array_temp__7_92_1 = SSAValue16;
					SSAValue8 = (parallel_ir_array_temp__7_92_1 == parfor_index_2_82);
					parallel_ir_array_temp__21_94_1 = SSAValue8;
					if (parallel_ir_array_temp__21_94_1)
					{
						;
						parallel_ir_array_temp__9_97_1 = _pa_rand_gen_arr.ARRAYELEM(parfor_index_1_82,((parfor_index_1_25) - (__hpat_loop_start_25)) + (1));
						HPAT__99_82.ARRAYELEM(parfor_index_1_82, parfor_index_2_82) = (HPAT__99_82.ARRAYELEM(parfor_index_1_82,parfor_index_2_82)+parallel_ir_array_temp__9_97_1);
					}
					;
					parallel_ir_array_temp__7_100_1 = SSAValue16;
					SSAValue9 = (parallel_ir_array_temp__7_100_1 == parfor_index_2_82);
					SSAValue26 = (SSAValue9) ? (1) : (0);
					HPAT__105_82.ARRAYELEM(parfor_index_1_82, parfor_index_2_82) = (HPAT__105_82.ARRAYELEM(parfor_index_1_82,parfor_index_2_82)+SSAValue26);
				}
			}
			;
		}
		;
		__hpat_reduce_95 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
		MPI_Allreduce(HPAT__99_82.data, __hpat_reduce_95.data, (10) * (numCenter), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
		HPAT__99_82 = __hpat_reduce_95;
		__hpat_reduce_96 = j2c_array<int64_t>::new_j2c_array_2d(NULL, 10, numCenter);
		MPI_Allreduce(HPAT__105_82.data, __hpat_reduce_96.data, (10) * (numCenter), MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);;
		HPAT__105_82 = __hpat_reduce_96;
		for(int64_t parfor_index_2_82 = 1; parfor_index_2_82 <= numCenter; parfor_index_2_82 += 1)
		{
			;
			for(int64_t parfor_index_1_82 = 1; parfor_index_1_82 <= 10; parfor_index_1_82 += 1)
			{
				;
				HPAT_SSAValue28_82.ARRAYELEM(parfor_index_1_82, parfor_index_2_82) = (double)HPAT__105_82.ARRAYELEM(parfor_index_1_82,parfor_index_2_82);
				HPAT_SSAValue2_82.ARRAYELEM(parfor_index_1_82, parfor_index_2_82) = (HPAT__99_82.ARRAYELEM(parfor_index_1_82,parfor_index_2_82)) / (HPAT_SSAValue28_82.ARRAYELEM(parfor_index_1_82,parfor_index_2_82));
				SSAValue4.ARRAYELEM(parfor_index_1_82,parfor_index_2_82) = HPAT_SSAValue2_82.ARRAYELEM(parfor_index_1_82,parfor_index_2_82);
			}
			;
		}
		;
		centroids = SSAValue4;
	}
	;
	0; if(__hpat_node_id==0) printf("exec time %lf\n", MPI_Wtime()-__hpat_t1);;
	*ret0 = centroids;
	return;

}


extern "C" void _ppkmeansp271_(int run_where, int64_t numCenter, int64_t iterNum, int64_t N ,  j2c_array< double > ** __restrict ret0 , bool genMain = true)
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
		mainFileData << numCenter << std::endl;
		mainFileData << iterNum << std::endl;
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
		mainFile << "    int64_t numCenter;" << std::endl;
		mainFile << "    mainFileData >> numCenter;" << std::endl;
		mainFile << "    int64_t iterNum;" << std::endl;
		mainFile << "    mainFileData >> iterNum;" << std::endl;
		mainFile << "    int64_t N;" << std::endl;
		mainFile << "    mainFileData >> N;" << std::endl;
		mainFile << "    _ppkmeansp271_(runwhere, numCenter, iterNum, N, &ret0, false);" << std::endl;
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

	ppkmeansp271(numCenter, iterNum, N, *ret0);
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
