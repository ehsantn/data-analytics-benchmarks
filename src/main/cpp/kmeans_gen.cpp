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
	j2c_array< double >  pp_pa_rand_gen_arr_17p283;
	int64_t SSAValue29;
	bool SSAValue33;
	int64_t SSAValue34;
	int64_t SSAValue36;
	bool SSAValue37;
	bool SSAValue13;
	j2c_array< int64_t >  SSAValue32;
	j2c_array< double >  SSAValue0;
	int64_t parfor_index_1_17;
	int64_t parfor_index_2_17;
	double SSAValue1;
	int64_t parfor_index_1_21;
	int64_t parfor_index_2_21;
	double SSAValue2;
	int64_t parfor_index_1_25;
	int64_t ppptemppp299;
	int64_t ppptempp_5p287;
	int64_t ppptempp_6p288;
	int64_t ppjp300;
	int64_t SSAValue3;
	j2c_array< double >  SSAValue4;
	double SSAValue5;
	int64_t SSAValue6;
	int64_t SSAValue7;
	bool SSAValue8;
	int64_t SSAValue9;
	int64_t SSAValue10;
	bool SSAValue11;
	bool SSAValue12;
	double parallel_ir_temp__13_parallel_ir_range__41Int64_skip__43Int64_1_1_1;
	double parallel_ir_temp__11_parallel_ir_range__47Int64_skip__49Int64_1_1_1;
	int64_t parfor_index_1_46;
	double SSAValue15;
	j2c_array< double >  parallel_ir_new_array_name_46_1;
	double parallel_ir_array_temp_SSAValue30_58_1;
	int64_t parfor_index_1_57;
	int64_t parallel_ir_save_array_len_1_57;
	int32_t SSAValue16;
	double SSAValue17;
	double parallel_ir_array_temp_SSAValue31_63_1;
	double parallel_ir_reduction_output_61;
	ParallelAcceleratorAPILibMTFloat64 m;
	int64_t SSAValue19;
	int64_t SSAValue20;
	int64_t ppip301;
	ParallelAcceleratorAPILibMTFloat64 parallel_ir_reduction_input_70_1;
	ParallelAcceleratorAPILibMTFloat64 pptemp_neutral_valp302;
	double SSAValue21;
	int64_t SSAValue22;
	double SSAValue23;
	double SSAValue24;
	bool SSAValue25;
	double SSAValue26;
	double SSAValue27;
	double SSAValue28;
	bool SSAValue30;
	int64_t parfor_index_1_75;
	int64_t parfor_index_2_75;
	double SSAValue35;
	double SSAValue38;
	int64_t parallel_ir_array_temp__7_85_1;
	int64_t parfor_index_1_84;
	bool SSAValue14;
	bool parallel_ir_array_temp__21_87_1;
	double parallel_ir_array_temp__9_90_1;
	double parallel_ir_reduction_output_88;
	int64_t parallel_ir_array_temp__7_93_1;
	bool SSAValue18;
	int64_t SSAValue31;
	int64_t parallel_ir_array_temp_SSAValue6_102_1;
	int64_t parallel_ir_reduction_output_100;
	j2c_array< double >  HPAT__103_75;
	j2c_array< int64_t >  HPAT__109_75;
	j2c_array< double >  HPAT_SSAValue35_75;
	j2c_array< double >  HPAT_SSAValue38_75;
	int32_t __hpat_num_pes;
	int32_t __hpat_node_id;
	int64_t __hpat_loop_start_17;
	int64_t __hpat_loop_end_17;
	int64_t __hpat_loop_div_17;
	int64_t __hpat_bcast_size_107;
	int64_t __hpat_dist_arr_start_87;
	int64_t __hpat_dist_arr_div_87;
	int64_t __hpat_dist_arr_count_87;
	int64_t __hpat_loop_start_25;
	int64_t __hpat_loop_end_25;
	int64_t __hpat_loop_div_25;
	j2c_array< double >  __hpat_reduce_88;
	j2c_array< int64_t >  __hpat_reduce_89;
	std::random_device cgen_rand_device;
	std::uniform_real_distribution<double> cgen_distribution(0.0,1.0);
	std::normal_distribution<double> cgen_n_distribution(0.0,1.0);
	std::default_random_engine cgen_rand_generator(cgen_rand_device());
	;;
	MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
	MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
	SSAValue33 = (1) <= (iterNum);
	SSAValue34 = (1) - (1);
	SSAValue29 = (SSAValue33) ? (iterNum) : (SSAValue34);
	SSAValue36 = (SSAValue29) + (1);
	__hpat_loop_div_17 = (N) / (__hpat_num_pes);
	__hpat_loop_start_17 = ((__hpat_node_id) * (__hpat_loop_div_17)) + (1);
	__hpat_loop_end_17 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_17);
	_pa_rand_gen_arr = j2c_array<double>::new_j2c_array_2d(NULL, 10, N);
	#pragma simd
	for ( parfor_index_2_17 = __hpat_loop_start_17; parfor_index_2_17 <= (int64_t)__hpat_loop_end_17; parfor_index_2_17 += 1)
	{
		for ( parfor_index_1_17 = 1; parfor_index_1_17 <= (int64_t)10; parfor_index_1_17 += 1)
		{
			;
			SSAValue1 = cgen_distribution(cgen_rand_generator);
			;
			_pa_rand_gen_arr.ARRAYELEM(parfor_index_1_17,((parfor_index_2_17) - (__hpat_loop_start_17)) + (1)) = SSAValue1;
		}
	}
	;
	if ((__hpat_node_id == 0))
	{
		;
		pp_pa_rand_gen_arr_17p283 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
		for ( parfor_index_2_21 = 1; parfor_index_2_21 <= (int64_t)numCenter; parfor_index_2_21 += 1)
		{
			for ( parfor_index_1_21 = 1; parfor_index_1_21 <= (int64_t)10; parfor_index_1_21 += 1)
			{
				;
				SSAValue2 = cgen_distribution(cgen_rand_generator);
				;
				pp_pa_rand_gen_arr_17p283.ARRAYELEM(parfor_index_1_21,parfor_index_2_21) = SSAValue2;
			}
		}
		;
	}
	;
	__hpat_bcast_size_107 = (10) * (numCenter);
	MPI_Bcast(pp_pa_rand_gen_arr_17p283.data, __hpat_bcast_size_107, MPI_DOUBLE, 0, MPI_COMM_WORLD);;
	centroids = pp_pa_rand_gen_arr_17p283;
	0; double __hpat_t1 = MPI_Wtime();
	ptempp = 1;
	while (1)
	{
		;
		SSAValue37 = (ptempp == SSAValue36);
		SSAValue13 = !(SSAValue37);
		if (!(SSAValue13)) break;
		ptempp = (ptempp) + (1);
		__hpat_dist_arr_div_87 = (N) / (__hpat_num_pes);
		__hpat_dist_arr_start_87 = (__hpat_node_id) * (__hpat_dist_arr_div_87);
		__hpat_dist_arr_count_87 = ((__hpat_node_id==__hpat_num_pes-1) ? N-__hpat_node_id*__hpat_dist_arr_div_87 : __hpat_dist_arr_div_87);
		SSAValue32 = j2c_array<int64_t>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_87);
		SSAValue0 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
		__hpat_loop_div_25 = (N) / (__hpat_num_pes);
		__hpat_loop_start_25 = ((__hpat_node_id) * (__hpat_loop_div_25)) + (1);
		__hpat_loop_end_25 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_25);
		HPAT__103_75 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
		HPAT__109_75 = j2c_array<int64_t>::new_j2c_array_2d(NULL, 10, numCenter);
		HPAT_SSAValue35_75 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
		HPAT_SSAValue38_75 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
		for(int64_t parfor_index_2_75 = 1; parfor_index_2_75 <= numCenter; parfor_index_2_75 += 1)
		{
			;
			for(int64_t parfor_index_1_75 = 1; parfor_index_1_75 <= 10; parfor_index_1_75 += 1)
			{
				;
				HPAT__103_75.ARRAYELEM(parfor_index_1_75, parfor_index_2_75) = 0.0;
				HPAT__109_75.ARRAYELEM(parfor_index_1_75, parfor_index_2_75) = 0;
			}
			;
		}
		;
		#pragma simd
		for ( parfor_index_1_25 = __hpat_loop_start_25; parfor_index_1_25 <= (int64_t)__hpat_loop_end_25; parfor_index_1_25 += 1)
		{
			;
			SSAValue8 = (1) <= (numCenter);
			SSAValue9 = (1) - (1);
			SSAValue6 = (SSAValue8) ? (numCenter) : (SSAValue9);
			SSAValue10 = (SSAValue6) - (1);
			SSAValue3 = (SSAValue10) + (1);
			SSAValue4 = j2c_array<double>::new_j2c_array_1d(NULL, ((1) <= (numCenter)) ? (numCenter) : (0));
			ppptemppp299 = 1;
			ppptempp_5p287 = 1;
			ppptempp_6p288 = 0;
			while (1)
			{
				;
				SSAValue11 = (ppptempp_6p288 == SSAValue3);
				SSAValue12 = !(SSAValue11);
				if (!(SSAValue12)) break;
				ppptempp_6p288 = (ppptempp_6p288) + (1);
				SSAValue7 = (ppptempp_5p287) + (1);
				ppjp300 = ppptempp_5p287;
				ppptempp_5p287 = SSAValue7;
				parallel_ir_new_array_name_46_1 = j2c_array<double>::new_j2c_array_1d(NULL, 10);
				for ( parfor_index_1_46 = 1; parfor_index_1_46 <= (int64_t)10; parfor_index_1_46 += 1)
				{
					;
					parallel_ir_temp__13_parallel_ir_range__41Int64_skip__43Int64_1_1_1 = _pa_rand_gen_arr.ARRAYELEM(parfor_index_1_46,((parfor_index_1_25) - (__hpat_loop_start_25)) + (1));
					parallel_ir_temp__11_parallel_ir_range__47Int64_skip__49Int64_1_1_1 = centroids.ARRAYELEM(parfor_index_1_46,ppjp300);
					SSAValue15 = (parallel_ir_temp__13_parallel_ir_range__41Int64_skip__43Int64_1_1_1) - (parallel_ir_temp__11_parallel_ir_range__47Int64_skip__49Int64_1_1_1);
					parallel_ir_new_array_name_46_1.ARRAYELEM(parfor_index_1_46) = SSAValue15;
				}
				;
				parallel_ir_save_array_len_1_57 = 10;
				parallel_ir_reduction_output_61 = 0.0;
				for ( parfor_index_1_57 = 1; parfor_index_1_57 <= (int64_t)parallel_ir_save_array_len_1_57; parfor_index_1_57 += 1)
				{
					;
					parallel_ir_array_temp_SSAValue30_58_1 = parallel_ir_new_array_name_46_1.ARRAYELEM(parfor_index_1_57);
					SSAValue16 = (int32_t)(2);
					SSAValue17 = pow(parallel_ir_array_temp_SSAValue30_58_1, SSAValue16);
					parallel_ir_reduction_output_61 = (parallel_ir_reduction_output_61+SSAValue17);
				}
				;
				SSAValue5 = sqrt(parallel_ir_reduction_output_61);;
				SSAValue4.ARRAYELEM(ppptemppp299) = SSAValue5;
				ppptemppp299 = (ppptemppp299) + (1);
			}
			;
			m = ParallelAcceleratorAPILibMTFloat64{SSAValue4.ARRAYELEM(1), 1};
			SSAValue19 = SSAValue4.ARRAYLEN();
			for ( ppip301 = 1; ppip301 <= (int64_t)SSAValue19; ppip301 += 1)
			{
				;
				SSAValue27 = SSAValue4.ARRAYELEM(ppip301);
				SSAValue28 = m.val;
				SSAValue30 = (SSAValue27) < (SSAValue28);
				if (SSAValue30)
				{
					;
					SSAValue26 = SSAValue4.ARRAYELEM(ppip301);
					m.val = SSAValue26;
					m.idx = ppip301;
				}
				;
			}
			;
			SSAValue20 = m.idx;
			for ( parfor_index_2_75 = 1; parfor_index_2_75 <= (int64_t)numCenter; parfor_index_2_75 += 1)
			{
				for ( parfor_index_1_75 = 1; parfor_index_1_75 <= (int64_t)10; parfor_index_1_75 += 1)
				{
					;
					parallel_ir_array_temp__7_85_1 = SSAValue20;
					SSAValue14 = (parallel_ir_array_temp__7_85_1 == parfor_index_2_75);
					parallel_ir_array_temp__21_87_1 = SSAValue14;
					if (parallel_ir_array_temp__21_87_1)
					{
						;
						parallel_ir_array_temp__9_90_1 = _pa_rand_gen_arr.ARRAYELEM(parfor_index_1_75,((parfor_index_1_25) - (__hpat_loop_start_25)) + (1));
						HPAT__103_75.ARRAYELEM(parfor_index_1_75, parfor_index_2_75) = (HPAT__103_75.ARRAYELEM(parfor_index_1_75,parfor_index_2_75)+parallel_ir_array_temp__9_90_1);
					}
					;
					parallel_ir_array_temp__7_93_1 = SSAValue20;
					SSAValue18 = (parallel_ir_array_temp__7_93_1 == parfor_index_2_75);
					SSAValue31 = (SSAValue18) ? (1) : (0);
					HPAT__109_75.ARRAYELEM(parfor_index_1_75, parfor_index_2_75) = (HPAT__109_75.ARRAYELEM(parfor_index_1_75,parfor_index_2_75)+SSAValue31);
				}
			}
			;
		}
		;
		__hpat_reduce_88 = j2c_array<double>::new_j2c_array_2d(NULL, 10, numCenter);
		MPI_Allreduce(HPAT__103_75.data, __hpat_reduce_88.data, (10) * (numCenter), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
		HPAT__103_75 = __hpat_reduce_88;
		__hpat_reduce_89 = j2c_array<int64_t>::new_j2c_array_2d(NULL, 10, numCenter);
		MPI_Allreduce(HPAT__109_75.data, __hpat_reduce_89.data, (10) * (numCenter), MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);;
		HPAT__109_75 = __hpat_reduce_89;
		for(int64_t parfor_index_2_75 = 1; parfor_index_2_75 <= numCenter; parfor_index_2_75 += 1)
		{
			;
			for(int64_t parfor_index_1_75 = 1; parfor_index_1_75 <= 10; parfor_index_1_75 += 1)
			{
				;
				HPAT_SSAValue35_75.ARRAYELEM(parfor_index_1_75, parfor_index_2_75) = (double)HPAT__109_75.ARRAYELEM(parfor_index_1_75,parfor_index_2_75);
				HPAT_SSAValue38_75.ARRAYELEM(parfor_index_1_75, parfor_index_2_75) = (HPAT__103_75.ARRAYELEM(parfor_index_1_75,parfor_index_2_75)) / (HPAT_SSAValue35_75.ARRAYELEM(parfor_index_1_75,parfor_index_2_75));
				SSAValue0.ARRAYELEM(parfor_index_1_75,parfor_index_2_75) = HPAT_SSAValue38_75.ARRAYELEM(parfor_index_1_75,parfor_index_2_75);
			}
			;
		}
		;
		centroids = SSAValue0;
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
