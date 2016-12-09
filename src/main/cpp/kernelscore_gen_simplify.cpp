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
j2c_array< double >  _Base_vect(double X1, double X2, double X3);
j2c_array< double >  _Base_vect(double X1, double X2, double X3)
{
	TupleFloat64Float64Float64 X = {X1, X2, X3};
	Basepvect pselfp;
	int64_t ptempp;
	int64_t ppptempp_4p293;
	int64_t ppptempp_5p294;
	int64_t i;
	int64_t ppptempp_7p295;
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
	ppptempp_4p293 = 1;
	ppptempp_5p294 = 0;
	label8 : ;
	if (!(!((ppptempp_5p294 == SSAValue1)))) goto label28;
	SSAValue3 = (ppptempp_5p294) + (1);
	ppptempp_5p294 = SSAValue3;
	SSAValue11 = ppptempp_4p293;
	SSAValue12 = (ppptempp_4p293) + (1);
	ppptempp_7p295 = 1;
	SSAValue13 = (1) + (1);
	i = SSAValue11;
	ppptempp_7p295 = SSAValue13;
	SSAValue14 = (2) + (1);
	ppptempp_4p293 = SSAValue12;
	ppptempp_7p295 = SSAValue14;
	SSAValue7 = ((double *)&X)[i - 1];
	SSAValue2.ARRAYELEM(ptempp) = SSAValue7;
	ptempp = (ptempp) + (1);
	label26 : ;
	goto label8;
	label28 : ;
	return SSAValue2;

}


void ppkernelscore_testp271(int64_t n, double * __restrict ret0)
{
	pppkernelscore_testp271 pselfp;
	j2c_array< double >  X;
	j2c_array< double >  points;
	int64_t N;
	double b;
	double exps;
	j2c_array< double >  arr;
	int64_t SSAValue11;
	bool SSAValue13;
	int64_t SSAValue14;
	int64_t SSAValue16;
	int64_t SSAValue18;
	int64_t parfor_index_1_3;
	int64_t parallel_ir_save_array_len_1_3;
	double SSAValue19;
	double parallel_ir_array_temp__10_6_2;
	int64_t ppip273p280;
	double parallel_ir_reduction_input_7_1;
	double pptemp_neutral_valp281;
	j2c_array< double >  d;
	double m;
	j2c_array< double >  SSAValue20;
	double SSAValue21;
	double SSAValue22;
	double SSAValue23;
	double SSAValue24;
	double SSAValue25;
	double SSAValue16pp1;
	j2c_array< double >  SSAValue26;
	j2c_array< double >  SSAValue27;
	int32_t SSAValue0;
	double SSAValue1;
	double SSAValue25pp2;
	double SSAValue2;
	double SSAValue3;
	j2c_array< double >  SSAValue4;
	j2c_array< double >  SSAValue5;
	double SSAValue6;
	double parallel_ir_array_temp__7_18_1;
	int64_t parfor_index_1_17;
	int64_t parallel_ir_save_array_len_1_17;
	double SSAValue7;
	j2c_array< double >  parallel_ir_new_array_name_17_1;
	double parallel_ir_array_temp__25_20_1;
	double parallel_ir_array_temp_SSAValue17_22_1;
	int64_t parfor_index_1_21;
	int64_t parallel_ir_save_array_len_1_21;
	int32_t SSAValue8;
	double SSAValue9;
	double parallel_ir_array_temp_SSAValue17_24_2;
	double parallel_ir_array_temp_SSAValue18_26_1;
	int64_t parfor_index_1_25;
	int64_t parallel_ir_save_array_len_1_25;
	double SSAValue10;
	double parallel_ir_array_temp_SSAValue18_28_2;
	double parallel_ir_array_temp_SSAValue2_30_1;
	int64_t parfor_index_1_29;
	int64_t parallel_ir_save_array_len_1_29;
	double SSAValue12;
	double parallel_ir_array_temp_SSAValue2_32_2;
	int64_t parfor_index_1_33;
	double parallel_ir_array_temp__4_35_1;
	int64_t parallel_ir_save_array_len_1_33;
	double parallel_ir_reduction_output_33;
	double parallel_ir_array_temp__4_38_1;
	int64_t parfor_index_1_37;
	int64_t parallel_ir_save_array_len_1_37;
	double SSAValue15;
	double parallel_ir_array_temp__4_40_2;
	double parallel_ir_array_temp_SSAValue31_42_1;
	int64_t parfor_index_1_41;
	int64_t parallel_ir_save_array_len_1_41;
	double SSAValue17;
	double parallel_ir_array_temp_SSAValue31_44_2;
	int64_t parfor_index_1_45;
	double parallel_ir_array_temp_SSAValue32_47_1;
	int64_t parallel_ir_save_array_len_1_45;
	double parallel_ir_reduction_output_45;
	int32_t __hpat_num_pes;
	int32_t __hpat_node_id;
	int64_t __hpat_dist_arr_start_1;
	int64_t __hpat_dist_arr_div_1;
	int64_t __hpat_dist_arr_count_1;
	int64_t __hpat_loop_start_3;
	int64_t __hpat_loop_end_3;
	int64_t __hpat_loop_div_3;
	int64_t __hpat_loop_start_7;
	int64_t __hpat_loop_end_7;
	int64_t __hpat_loop_div_7;
	double __hpat_reduce_2;
	std::random_device cgen_rand_device;
	std::uniform_real_distribution<double> cgen_distribution(0.0,1.0);
	std::normal_distribution<double> cgen_n_distribution(0.0,1.0);
	std::default_random_engine cgen_rand_generator(cgen_rand_device());
	;;
	MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
	MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
	SSAValue14 = (1) - (1);
	__hpat_dist_arr_div_1 = (n) / (__hpat_num_pes);
	__hpat_dist_arr_start_1 = (__hpat_node_id) * (__hpat_dist_arr_div_1);
	__hpat_dist_arr_count_1 = ((__hpat_node_id==__hpat_num_pes-1) ? n-__hpat_node_id*__hpat_dist_arr_div_1 : __hpat_dist_arr_div_1);
	arr = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
	parallel_ir_save_array_len_1_3 = n;
	__hpat_loop_div_3 = (parallel_ir_save_array_len_1_3) / (__hpat_num_pes);
	__hpat_loop_start_3 = ((__hpat_node_id) * (__hpat_loop_div_3)) + (1);
	__hpat_loop_end_3 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_3 : (__hpat_node_id+1)*__hpat_loop_div_3);
	for ( parfor_index_1_3 = __hpat_loop_start_3; parfor_index_1_3 <= (int64_t)__hpat_loop_end_3; parfor_index_1_3 += 1)
	{
		;
		SSAValue19 = cgen_distribution(cgen_rand_generator);
		;
		parallel_ir_array_temp__10_6_2 = SSAValue19;
		arr.ARRAYELEM(((parfor_index_1_3) - (__hpat_loop_start_3)) + (1)) = parallel_ir_array_temp__10_6_2;
	}
	t1 = MPI_Wtime();
	points = _Base_vect(-1.0,2.0,5.0);
	exps = 0.0;
	SSAValue13 = (1) <= (n);
	SSAValue11 = (SSAValue13) ? (n) : (SSAValue14);
	SSAValue16 = (SSAValue11) - (1);
	SSAValue18 = (SSAValue16) + (1);
	__hpat_loop_div_7 = (SSAValue18) / (__hpat_num_pes);
	__hpat_loop_start_7 = ((__hpat_node_id) * (__hpat_loop_div_7)) + (1);
	__hpat_loop_end_7 = ((__hpat_node_id==__hpat_num_pes-1) ? SSAValue18 : (__hpat_node_id+1)*__hpat_loop_div_7);
//	#pragma simd reduction(+: exps)
	for ( ppip273p280 = __hpat_loop_start_7; ppip273p280 <= (int64_t)__hpat_loop_end_7; ppip273p280 += 1)
	{
		;
		SSAValue0 = (int32_t)(2);
		SSAValue21 = pow(0.5, SSAValue0);
		SSAValue1 = (double)2;
		SSAValue25pp2 = (SSAValue1) * (SSAValue21);
		SSAValue2 = (double)3;
		SSAValue3 = (0.5) * (SSAValue2);
		SSAValue22 = log(SSAValue3);;
		SSAValue16pp1 = arr.ARRAYELEM(((ppip273p280) - (__hpat_loop_start_7)) + (1));
		parallel_ir_save_array_len_1_17 = 3;
		parallel_ir_new_array_name_17_1 = j2c_array<double>::new_j2c_array_1d(NULL, 3);
		for ( parfor_index_1_17 = 1; parfor_index_1_17 <= (int64_t)parallel_ir_save_array_len_1_17; parfor_index_1_17 += 1)
		{
			;
			parallel_ir_array_temp__7_18_1 = points.ARRAYELEM(parfor_index_1_17);
			SSAValue7 = (SSAValue16pp1) - (parallel_ir_array_temp__7_18_1);
			parallel_ir_array_temp__25_20_1 = SSAValue7;
			parallel_ir_new_array_name_17_1.ARRAYELEM(parfor_index_1_17) = parallel_ir_array_temp__25_20_1;
		}
		;
		SSAValue26 = parallel_ir_new_array_name_17_1;
		parallel_ir_save_array_len_1_21 = 3;
		for ( parfor_index_1_21 = 1; parfor_index_1_21 <= (int64_t)parallel_ir_save_array_len_1_21; parfor_index_1_21 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue17_22_1 = SSAValue26.ARRAYELEM(parfor_index_1_21);
			SSAValue8 = (int32_t)(2);
			SSAValue9 = pow(parallel_ir_array_temp_SSAValue17_22_1, SSAValue8);
			parallel_ir_array_temp_SSAValue17_24_2 = SSAValue9;
			SSAValue26.ARRAYELEM(parfor_index_1_21) = parallel_ir_array_temp_SSAValue17_24_2;
		}
		;
		SSAValue27 = SSAValue26;
		parallel_ir_save_array_len_1_25 = 3;
		for ( parfor_index_1_25 = 1; parfor_index_1_25 <= (int64_t)parallel_ir_save_array_len_1_25; parfor_index_1_25 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue18_26_1 = SSAValue27.ARRAYELEM(parfor_index_1_25);
			SSAValue10 = -(parallel_ir_array_temp_SSAValue18_26_1);
			parallel_ir_array_temp_SSAValue18_28_2 = SSAValue10;
			SSAValue27.ARRAYELEM(parfor_index_1_25) = parallel_ir_array_temp_SSAValue18_28_2;
		}
		;
		SSAValue20 = SSAValue27;
		parallel_ir_save_array_len_1_29 = 3;
		for ( parfor_index_1_29 = 1; parfor_index_1_29 <= (int64_t)parallel_ir_save_array_len_1_29; parfor_index_1_29 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue2_30_1 = SSAValue20.ARRAYELEM(parfor_index_1_29);
			SSAValue12 = (parallel_ir_array_temp_SSAValue2_30_1) / (SSAValue25pp2);
			parallel_ir_array_temp_SSAValue2_32_2 = SSAValue12;
			SSAValue20.ARRAYELEM(parfor_index_1_29) = parallel_ir_array_temp_SSAValue2_32_2;
		}
		;
		d = SSAValue20;
		parallel_ir_save_array_len_1_33 = d.ARRAYSIZE(1);
		parallel_ir_reduction_output_33 = DBL_MAX;
		for ( parfor_index_1_33 = 1; parfor_index_1_33 <= (int64_t)parallel_ir_save_array_len_1_33; parfor_index_1_33 += 1)
		{
			;
			parallel_ir_array_temp__4_35_1 = d.ARRAYELEM(parfor_index_1_33);
			parallel_ir_reduction_output_33 = std::min((double)parallel_ir_reduction_output_33,(double)parallel_ir_array_temp__4_35_1);
		}
		;
		m = parallel_ir_reduction_output_33;
		SSAValue23 = (m) - (SSAValue22);
		parallel_ir_save_array_len_1_37 = 3;
		for ( parfor_index_1_37 = 1; parfor_index_1_37 <= (int64_t)parallel_ir_save_array_len_1_37; parfor_index_1_37 += 1)
		{
			;
			parallel_ir_array_temp__4_38_1 = d.ARRAYELEM(parfor_index_1_37);
			SSAValue15 = (parallel_ir_array_temp__4_38_1) - (m);
			parallel_ir_array_temp__4_40_2 = SSAValue15;
			d.ARRAYELEM(parfor_index_1_37) = parallel_ir_array_temp__4_40_2;
		}
		;
		SSAValue4 = d;
		parallel_ir_save_array_len_1_41 = 3;
		for ( parfor_index_1_41 = 1; parfor_index_1_41 <= (int64_t)parallel_ir_save_array_len_1_41; parfor_index_1_41 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue31_42_1 = SSAValue4.ARRAYELEM(parfor_index_1_41);
			SSAValue17 = exp(parallel_ir_array_temp_SSAValue31_42_1);
			parallel_ir_array_temp_SSAValue31_44_2 = SSAValue17;
			SSAValue4.ARRAYELEM(parfor_index_1_41) = parallel_ir_array_temp_SSAValue31_44_2;
		}
		;
		SSAValue5 = SSAValue4;
		parallel_ir_save_array_len_1_45 = SSAValue5.ARRAYSIZE(1);
		parallel_ir_reduction_output_45 = 0.0;
		for ( parfor_index_1_45 = 1; parfor_index_1_45 <= (int64_t)parallel_ir_save_array_len_1_45; parfor_index_1_45 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue32_47_1 = SSAValue5.ARRAYELEM(parfor_index_1_45);
			parallel_ir_reduction_output_45 = (parallel_ir_reduction_output_45+parallel_ir_array_temp_SSAValue32_47_1);
		}
		;
		SSAValue6 = parallel_ir_reduction_output_45;
		SSAValue24 = log(SSAValue6);;
		SSAValue25 = (SSAValue23) + (SSAValue24);
		exps = (exps) + (SSAValue25);
	}
	;
	__hpat_reduce_2 = 0;
	MPI_Allreduce(&exps, &__hpat_reduce_2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
	exps = __hpat_reduce_2;
	t2 = MPI_Wtime();
	*ret0 = exps;
	return;

}


void ppkernelscore_testp271_unaliased(int64_t n, double * __restrict ret0)
{
	pppkernelscore_testp271 pselfp;
	j2c_array< double >  X;
	j2c_array< double >  points;
	int64_t N;
	double b;
	double exps;
	j2c_array< double >  arr;
	int64_t SSAValue11;
	bool SSAValue13;
	int64_t SSAValue14;
	int64_t SSAValue16;
	int64_t SSAValue18;
	int64_t parfor_index_1_3;
	int64_t parallel_ir_save_array_len_1_3;
	double SSAValue19;
	double parallel_ir_array_temp__10_6_2;
	int64_t ppip273p280;
	double parallel_ir_reduction_input_7_1;
	double pptemp_neutral_valp281;
	j2c_array< double >  d;
	double m;
	j2c_array< double >  SSAValue20;
	double SSAValue21;
	double SSAValue22;
	double SSAValue23;
	double SSAValue24;
	double SSAValue25;
	double SSAValue16pp1;
	j2c_array< double >  SSAValue26;
	j2c_array< double >  SSAValue27;
	int32_t SSAValue0;
	double SSAValue1;
	double SSAValue25pp2;
	double SSAValue2;
	double SSAValue3;
	j2c_array< double >  SSAValue4;
	j2c_array< double >  SSAValue5;
	double SSAValue6;
	double parallel_ir_array_temp__7_18_1;
	int64_t parfor_index_1_17;
	int64_t parallel_ir_save_array_len_1_17;
	double SSAValue7;
	j2c_array< double >  parallel_ir_new_array_name_17_1;
	double parallel_ir_array_temp__25_20_1;
	double parallel_ir_array_temp_SSAValue17_22_1;
	int64_t parfor_index_1_21;
	int64_t parallel_ir_save_array_len_1_21;
	int32_t SSAValue8;
	double SSAValue9;
	double parallel_ir_array_temp_SSAValue17_24_2;
	double parallel_ir_array_temp_SSAValue18_26_1;
	int64_t parfor_index_1_25;
	int64_t parallel_ir_save_array_len_1_25;
	double SSAValue10;
	double parallel_ir_array_temp_SSAValue18_28_2;
	double parallel_ir_array_temp_SSAValue2_30_1;
	int64_t parfor_index_1_29;
	int64_t parallel_ir_save_array_len_1_29;
	double SSAValue12;
	double parallel_ir_array_temp_SSAValue2_32_2;
	int64_t parfor_index_1_33;
	double parallel_ir_array_temp__4_35_1;
	int64_t parallel_ir_save_array_len_1_33;
	double parallel_ir_reduction_output_33;
	double parallel_ir_array_temp__4_38_1;
	int64_t parfor_index_1_37;
	int64_t parallel_ir_save_array_len_1_37;
	double SSAValue15;
	double parallel_ir_array_temp__4_40_2;
	double parallel_ir_array_temp_SSAValue31_42_1;
	int64_t parfor_index_1_41;
	int64_t parallel_ir_save_array_len_1_41;
	double SSAValue17;
	double parallel_ir_array_temp_SSAValue31_44_2;
	int64_t parfor_index_1_45;
	double parallel_ir_array_temp_SSAValue32_47_1;
	int64_t parallel_ir_save_array_len_1_45;
	double parallel_ir_reduction_output_45;
	int32_t __hpat_num_pes;
	int32_t __hpat_node_id;
	int64_t __hpat_dist_arr_start_1;
	int64_t __hpat_dist_arr_div_1;
	int64_t __hpat_dist_arr_count_1;
	int64_t __hpat_loop_start_3;
	int64_t __hpat_loop_end_3;
	int64_t __hpat_loop_div_3;
	int64_t __hpat_loop_start_7;
	int64_t __hpat_loop_end_7;
	int64_t __hpat_loop_div_7;
	double __hpat_reduce_2;
	std::random_device cgen_rand_device;
	std::uniform_real_distribution<double> cgen_distribution(0.0,1.0);
	std::normal_distribution<double> cgen_n_distribution(0.0,1.0);
	std::default_random_engine cgen_rand_generator(cgen_rand_device());
	;;
	MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
	MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
	SSAValue14 = (1) - (1);
	__hpat_dist_arr_div_1 = (n) / (__hpat_num_pes);
	__hpat_dist_arr_start_1 = (__hpat_node_id) * (__hpat_dist_arr_div_1);
	__hpat_dist_arr_count_1 = ((__hpat_node_id==__hpat_num_pes-1) ? n-__hpat_node_id*__hpat_dist_arr_div_1 : __hpat_dist_arr_div_1);
	arr = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
	parallel_ir_save_array_len_1_3 = n;
	__hpat_loop_div_3 = (parallel_ir_save_array_len_1_3) / (__hpat_num_pes);
	__hpat_loop_start_3 = ((__hpat_node_id) * (__hpat_loop_div_3)) + (1);
	__hpat_loop_end_3 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_3 : (__hpat_node_id+1)*__hpat_loop_div_3);
	for ( parfor_index_1_3 = __hpat_loop_start_3; parfor_index_1_3 <= (int64_t)__hpat_loop_end_3; parfor_index_1_3 += 1)
	{
		;
		SSAValue19 = cgen_distribution(cgen_rand_generator);
		;
		parallel_ir_array_temp__10_6_2 = SSAValue19;
		arr.ARRAYELEM(((parfor_index_1_3) - (__hpat_loop_start_3)) + (1)) = parallel_ir_array_temp__10_6_2;
	}
	;
	points = _Base_vect(-1.0,2.0,5.0);
	exps = 0.0;
	SSAValue13 = (1) <= (n);
	SSAValue11 = (SSAValue13) ? (n) : (SSAValue14);
	SSAValue16 = (SSAValue11) - (1);
	SSAValue18 = (SSAValue16) + (1);
	__hpat_loop_div_7 = (SSAValue18) / (__hpat_num_pes);
	__hpat_loop_start_7 = ((__hpat_node_id) * (__hpat_loop_div_7)) + (1);
	__hpat_loop_end_7 = ((__hpat_node_id==__hpat_num_pes-1) ? SSAValue18 : (__hpat_node_id+1)*__hpat_loop_div_7);
	#pragma simd reduction(+: exps)
	for ( ppip273p280 = __hpat_loop_start_7; ppip273p280 <= (int64_t)__hpat_loop_end_7; ppip273p280 += 1)
	{
		;
		SSAValue0 = (int32_t)(2);
		SSAValue21 = pow(0.5, SSAValue0);
		SSAValue1 = (double)2;
		SSAValue25pp2 = (SSAValue1) * (SSAValue21);
		SSAValue2 = (double)3;
		SSAValue3 = (0.5) * (SSAValue2);
		SSAValue22 = log(SSAValue3);;
		SSAValue16pp1 = arr.ARRAYELEM(((ppip273p280) - (__hpat_loop_start_7)) + (1));
		parallel_ir_save_array_len_1_17 = 3;
		parallel_ir_new_array_name_17_1 = j2c_array<double>::new_j2c_array_1d(NULL, 3);
		for ( parfor_index_1_17 = 1; parfor_index_1_17 <= (int64_t)parallel_ir_save_array_len_1_17; parfor_index_1_17 += 1)
		{
			;
			parallel_ir_array_temp__7_18_1 = points.ARRAYELEM(parfor_index_1_17);
			SSAValue7 = (SSAValue16pp1) - (parallel_ir_array_temp__7_18_1);
			parallel_ir_array_temp__25_20_1 = SSAValue7;
			parallel_ir_new_array_name_17_1.ARRAYELEM(parfor_index_1_17) = parallel_ir_array_temp__25_20_1;
		}
		;
		SSAValue26 = parallel_ir_new_array_name_17_1;
		parallel_ir_save_array_len_1_21 = 3;
		for ( parfor_index_1_21 = 1; parfor_index_1_21 <= (int64_t)parallel_ir_save_array_len_1_21; parfor_index_1_21 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue17_22_1 = SSAValue26.ARRAYELEM(parfor_index_1_21);
			SSAValue8 = (int32_t)(2);
			SSAValue9 = pow(parallel_ir_array_temp_SSAValue17_22_1, SSAValue8);
			parallel_ir_array_temp_SSAValue17_24_2 = SSAValue9;
			SSAValue26.ARRAYELEM(parfor_index_1_21) = parallel_ir_array_temp_SSAValue17_24_2;
		}
		;
		SSAValue27 = SSAValue26;
		parallel_ir_save_array_len_1_25 = 3;
		for ( parfor_index_1_25 = 1; parfor_index_1_25 <= (int64_t)parallel_ir_save_array_len_1_25; parfor_index_1_25 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue18_26_1 = SSAValue27.ARRAYELEM(parfor_index_1_25);
			SSAValue10 = -(parallel_ir_array_temp_SSAValue18_26_1);
			parallel_ir_array_temp_SSAValue18_28_2 = SSAValue10;
			SSAValue27.ARRAYELEM(parfor_index_1_25) = parallel_ir_array_temp_SSAValue18_28_2;
		}
		;
		SSAValue20 = SSAValue27;
		parallel_ir_save_array_len_1_29 = 3;
		for ( parfor_index_1_29 = 1; parfor_index_1_29 <= (int64_t)parallel_ir_save_array_len_1_29; parfor_index_1_29 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue2_30_1 = SSAValue20.ARRAYELEM(parfor_index_1_29);
			SSAValue12 = (parallel_ir_array_temp_SSAValue2_30_1) / (SSAValue25pp2);
			parallel_ir_array_temp_SSAValue2_32_2 = SSAValue12;
			SSAValue20.ARRAYELEM(parfor_index_1_29) = parallel_ir_array_temp_SSAValue2_32_2;
		}
		;
		d = SSAValue20;
		parallel_ir_save_array_len_1_33 = d.ARRAYSIZE(1);
		parallel_ir_reduction_output_33 = DBL_MAX;
		for ( parfor_index_1_33 = 1; parfor_index_1_33 <= (int64_t)parallel_ir_save_array_len_1_33; parfor_index_1_33 += 1)
		{
			;
			parallel_ir_array_temp__4_35_1 = d.ARRAYELEM(parfor_index_1_33);
			parallel_ir_reduction_output_33 = std::min((double)parallel_ir_reduction_output_33,(double)parallel_ir_array_temp__4_35_1);
		}
		;
		m = parallel_ir_reduction_output_33;
		SSAValue23 = (m) - (SSAValue22);
		parallel_ir_save_array_len_1_37 = 3;
		for ( parfor_index_1_37 = 1; parfor_index_1_37 <= (int64_t)parallel_ir_save_array_len_1_37; parfor_index_1_37 += 1)
		{
			;
			parallel_ir_array_temp__4_38_1 = d.ARRAYELEM(parfor_index_1_37);
			SSAValue15 = (parallel_ir_array_temp__4_38_1) - (m);
			parallel_ir_array_temp__4_40_2 = SSAValue15;
			d.ARRAYELEM(parfor_index_1_37) = parallel_ir_array_temp__4_40_2;
		}
		;
		SSAValue4 = d;
		parallel_ir_save_array_len_1_41 = 3;
		for ( parfor_index_1_41 = 1; parfor_index_1_41 <= (int64_t)parallel_ir_save_array_len_1_41; parfor_index_1_41 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue31_42_1 = SSAValue4.ARRAYELEM(parfor_index_1_41);
			SSAValue17 = exp(parallel_ir_array_temp_SSAValue31_42_1);
			parallel_ir_array_temp_SSAValue31_44_2 = SSAValue17;
			SSAValue4.ARRAYELEM(parfor_index_1_41) = parallel_ir_array_temp_SSAValue31_44_2;
		}
		;
		SSAValue5 = SSAValue4;
		parallel_ir_save_array_len_1_45 = SSAValue5.ARRAYSIZE(1);
		parallel_ir_reduction_output_45 = 0.0;
		for ( parfor_index_1_45 = 1; parfor_index_1_45 <= (int64_t)parallel_ir_save_array_len_1_45; parfor_index_1_45 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue32_47_1 = SSAValue5.ARRAYELEM(parfor_index_1_45);
			parallel_ir_reduction_output_45 = (parallel_ir_reduction_output_45+parallel_ir_array_temp_SSAValue32_47_1);
		}
		;
		SSAValue6 = parallel_ir_reduction_output_45;
		SSAValue24 = log(SSAValue6);;
		SSAValue25 = (SSAValue23) + (SSAValue24);
		exps = (exps) + (SSAValue25);
	}
	;
	__hpat_reduce_2 = 0;
	MPI_Allreduce(&exps, &__hpat_reduce_2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
	exps = __hpat_reduce_2;
	*ret0 = exps;
	return;

}


extern "C" void _ppkernelscore_testp271_unaliased_(int run_where, int64_t n , double* __restrict ret0 , bool genMain = true)
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
		mainFileData << n << std::endl;
		mainFileData.close();
		std::ofstream mainFile(newMain.str());
		mainFile << "#include \"" << __FILE__ << "\"" << std::endl;
		mainFile << "int main(int argc, char *argv[]) {" << std::endl;
		mainFile << "    MPI_Init(&argc, &argv);" << std::endl;
		mainFile << "    std::ifstream mainFileData(\"" << newMainData.str() << "\", std::ios::in | std::ios::binary);" << std::endl;
		mainFile << "    int runwhere;" << std::endl;
		mainFile << "    mainFileData >> runwhere;" << std::endl;
		mainFile << "    double ret0;" << std::endl;
		mainFile << "    int64_t n;" << std::endl;
		mainFile << "    mainFileData >> n;" << std::endl;
		mainFile << "    _ppkernelscore_testp271_unaliased_(runwhere, n, &ret0, false);" << std::endl;
		mainFile << "    MPI_Finalize();" << std::endl;
		mainFile << "    return 0;" << std::endl;
		mainFile << "}" << std::endl;
		mainFile.close();
		std::ofstream mainFileSh(newMainSh.str());
		mainFileSh << "#!/bin/sh" << std::endl;
		mainFileSh << "mpiicpc -O3 -std=c++11 -I/usr/local/hdf5/include  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << "  -lm " << std::endl;
		mainFileSh.close();
	}

	ppkernelscore_testp271_unaliased(n, ret0);
}


extern "C" void _ppkernelscore_testp271_(int run_where, int64_t n , double* __restrict ret0 , bool genMain = true)
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
		mainFileData << n << std::endl;
		mainFileData.close();
		std::ofstream mainFile(newMain.str());
		mainFile << "#include \"" << __FILE__ << "\"" << std::endl;
		mainFile << "int main(int argc, char *argv[]) {" << std::endl;
		mainFile << "    MPI_Init(&argc, &argv);" << std::endl;
		mainFile << "    std::ifstream mainFileData(\"" << newMainData.str() << "\", std::ios::in | std::ios::binary);" << std::endl;
		mainFile << "    int runwhere;" << std::endl;
		mainFile << "    mainFileData >> runwhere;" << std::endl;
		mainFile << "    double ret0;" << std::endl;
		mainFile << "    int64_t n;" << std::endl;
		mainFile << "    mainFileData >> n;" << std::endl;
		mainFile << "    _ppkernelscore_testp271_(runwhere, n, &ret0, false);" << std::endl;
		mainFile << "    MPI_Finalize();" << std::endl;
		mainFile << "    return 0;" << std::endl;
		mainFile << "}" << std::endl;
		mainFile.close();
		std::ofstream mainFileSh(newMainSh.str());
		mainFileSh << "#!/bin/sh" << std::endl;
		mainFileSh << "mpiicpc -O3 -std=c++11 -I/usr/local/hdf5/include  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << "  -lm " << std::endl;
		mainFileSh.close();
	}

	ppkernelscore_testp271(n, ret0);
}


extern "C"
void *j2c_array_new(int key, void*data, unsigned ndim, int64_t *dims)
{
	void *a = NULL;
	switch(key)
	{
		default:
			fprintf(stderr, "j2c_array_new called with invalid key %d", key);
			assert(false);
			break;
	}
	return a;
}
