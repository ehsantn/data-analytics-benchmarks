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
	int64_t SSAValue9;
	int64_t SSAValue11;
	bool SSAValue13;
	int64_t SSAValue14;
	int64_t SSAValue16;
	int64_t SSAValue18;
	int64_t parfor_index_1_1;
	int64_t parallel_ir_save_array_len_1_1;
	double SSAValue19;
	double parallel_ir_array_temp__10_4_2;
	int64_t ppip273p280;
	double parallel_ir_reduction_input_5_1;
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
	int32_t SSAValue28;
	double SSAValue0;
	double SSAValue25pp2;
	double SSAValue1;
	double SSAValue2;
	j2c_array< double >  SSAValue3;
	j2c_array< double >  SSAValue4;
	double SSAValue5;
	double parallel_ir_array_temp__7_16_1;
	int64_t parfor_index_1_15;
	int64_t parallel_ir_save_array_len_1_15;
	double SSAValue6;
	j2c_array< double >  parallel_ir_new_array_name_15_1;
	double parallel_ir_array_temp__23_18_1;
	double parallel_ir_array_temp_SSAValue17_20_1;
	int64_t parfor_index_1_19;
	int64_t parallel_ir_save_array_len_1_19;
	int32_t SSAValue7;
	double SSAValue8;
	double parallel_ir_array_temp_SSAValue17_22_2;
	double parallel_ir_array_temp_SSAValue18_24_1;
	int64_t parfor_index_1_23;
	int64_t parallel_ir_save_array_len_1_23;
	double SSAValue10;
	double parallel_ir_array_temp_SSAValue18_26_2;
	double parallel_ir_array_temp_SSAValue2_28_1;
	int64_t parfor_index_1_27;
	int64_t parallel_ir_save_array_len_1_27;
	double SSAValue12;
	double parallel_ir_array_temp_SSAValue2_30_2;
	int64_t parfor_index_1_31;
	double parallel_ir_array_temp__4_33_1;
	int64_t parallel_ir_save_array_len_1_31;
	double parallel_ir_reduction_output_31;
	double parallel_ir_array_temp__4_36_1;
	int64_t parfor_index_1_35;
	int64_t parallel_ir_save_array_len_1_35;
	double SSAValue15;
	double parallel_ir_array_temp__4_38_2;
	double parallel_ir_array_temp_SSAValue31_40_1;
	int64_t parfor_index_1_39;
	int64_t parallel_ir_save_array_len_1_39;
	double SSAValue17;
	double parallel_ir_array_temp_SSAValue31_42_2;
	int64_t parfor_index_1_43;
	double parallel_ir_array_temp_SSAValue32_45_1;
	int64_t parallel_ir_save_array_len_1_43;
	double parallel_ir_reduction_output_43;
	int32_t __hpat_num_pes;
	int32_t __hpat_node_id;
	int64_t __hpat_dist_arr_start_1;
	int64_t __hpat_dist_arr_div_1;
	int64_t __hpat_dist_arr_count_1;
	int64_t __hpat_loop_start_1;
	int64_t __hpat_loop_end_1;
	int64_t __hpat_loop_div_1;
	int64_t __hpat_loop_start_5;
	int64_t __hpat_loop_end_5;
	int64_t __hpat_loop_div_5;
	double __hpat_reduce_2;
	std::random_device cgen_rand_device;
	std::uniform_real_distribution<double> cgen_distribution(0.0,1.0);
	std::normal_distribution<double> cgen_n_distribution(0.0,1.0);
	std::default_random_engine cgen_rand_generator(cgen_rand_device());
	;;
	MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
	MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
	SSAValue14 = (1) - (1);
	exps = 0.0;
	__hpat_dist_arr_div_1 = (n) / (__hpat_num_pes);
	__hpat_dist_arr_start_1 = (__hpat_node_id) * (__hpat_dist_arr_div_1);
	__hpat_dist_arr_count_1 = ((__hpat_node_id==__hpat_num_pes-1) ? n-__hpat_node_id*__hpat_dist_arr_div_1 : __hpat_dist_arr_div_1);
	arr = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
	parallel_ir_save_array_len_1_1 = n;
	__hpat_loop_div_1 = (parallel_ir_save_array_len_1_1) / (__hpat_num_pes);
	__hpat_loop_start_1 = ((__hpat_node_id) * (__hpat_loop_div_1)) + (1);
	__hpat_loop_end_1 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_1 : (__hpat_node_id+1)*__hpat_loop_div_1);
	for ( parfor_index_1_1 = __hpat_loop_start_1; parfor_index_1_1 <= (int64_t)__hpat_loop_end_1; parfor_index_1_1 += 1)
	{
		;
		SSAValue19 = cgen_distribution(cgen_rand_generator);
		;
		parallel_ir_array_temp__10_4_2 = SSAValue19;
		arr.ARRAYELEM(((parfor_index_1_1) - (__hpat_loop_start_1)) + (1)) = parallel_ir_array_temp__10_4_2;
	}
	;
	X = arr;
    t1 = MPI_Wtime();
	points = _Base_vect(-1.0,2.0,5.0);
	N = points.ARRAYSIZE(1);
	b = 0.5;
	SSAValue9 = n;
	SSAValue13 = (1) <= (SSAValue9);
	SSAValue11 = (SSAValue13) ? (SSAValue9) : (SSAValue14);
	SSAValue16 = (SSAValue11) - (1);
	SSAValue18 = (SSAValue16) + (1);
	__hpat_loop_div_5 = (SSAValue18) / (__hpat_num_pes);
	__hpat_loop_start_5 = ((__hpat_node_id) * (__hpat_loop_div_5)) + (1);
	__hpat_loop_end_5 = ((__hpat_node_id==__hpat_num_pes-1) ? SSAValue18 : (__hpat_node_id+1)*__hpat_loop_div_5);
	//#pragma simd reduction(+: exps)
	for ( ppip273p280 = __hpat_loop_start_5; ppip273p280 <= (int64_t)__hpat_loop_end_5; ppip273p280 += 1)
	{
		;
		SSAValue28 = (int32_t)(2);
		SSAValue21 = pow(b, SSAValue28);
		SSAValue0 = (double)2;
		SSAValue25pp2 = (SSAValue0) * (SSAValue21);
		SSAValue1 = (double)N;
		SSAValue2 = (b) * (SSAValue1);
		SSAValue22 = log(SSAValue2);;
		SSAValue16pp1 = X.ARRAYELEM(((ppip273p280) - (__hpat_loop_start_5)) + (1));
		parallel_ir_save_array_len_1_15 = points.ARRAYSIZE(1);
		parallel_ir_new_array_name_15_1 = j2c_array<double>::new_j2c_array_1d(NULL, parallel_ir_save_array_len_1_15);
		for ( parfor_index_1_15 = 1; parfor_index_1_15 <= (int64_t)parallel_ir_save_array_len_1_15; parfor_index_1_15 += 1)
		{
			;
			parallel_ir_array_temp__7_16_1 = points.ARRAYELEM(parfor_index_1_15);
			SSAValue6 = (SSAValue16pp1) - (parallel_ir_array_temp__7_16_1);
			parallel_ir_array_temp__23_18_1 = SSAValue6;
			parallel_ir_new_array_name_15_1.ARRAYELEM(parfor_index_1_15) = parallel_ir_array_temp__23_18_1;
		}
		;
		SSAValue26 = parallel_ir_new_array_name_15_1;
		parallel_ir_save_array_len_1_19 = SSAValue26.ARRAYSIZE(1);
		for ( parfor_index_1_19 = 1; parfor_index_1_19 <= (int64_t)parallel_ir_save_array_len_1_19; parfor_index_1_19 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue17_20_1 = SSAValue26.ARRAYELEM(parfor_index_1_19);
			SSAValue7 = (int32_t)(2);
			SSAValue8 = pow(parallel_ir_array_temp_SSAValue17_20_1, SSAValue7);
			parallel_ir_array_temp_SSAValue17_22_2 = SSAValue8;
			SSAValue26.ARRAYELEM(parfor_index_1_19) = parallel_ir_array_temp_SSAValue17_22_2;
		}
		;
		SSAValue27 = SSAValue26;
		parallel_ir_save_array_len_1_23 = SSAValue27.ARRAYSIZE(1);
		for ( parfor_index_1_23 = 1; parfor_index_1_23 <= (int64_t)parallel_ir_save_array_len_1_23; parfor_index_1_23 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue18_24_1 = SSAValue27.ARRAYELEM(parfor_index_1_23);
			SSAValue10 = -(parallel_ir_array_temp_SSAValue18_24_1);
			parallel_ir_array_temp_SSAValue18_26_2 = SSAValue10;
			SSAValue27.ARRAYELEM(parfor_index_1_23) = parallel_ir_array_temp_SSAValue18_26_2;
		}
		;
		SSAValue20 = SSAValue27;
		parallel_ir_save_array_len_1_27 = SSAValue20.ARRAYSIZE(1);
		for ( parfor_index_1_27 = 1; parfor_index_1_27 <= (int64_t)parallel_ir_save_array_len_1_27; parfor_index_1_27 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue2_28_1 = SSAValue20.ARRAYELEM(parfor_index_1_27);
			SSAValue12 = (parallel_ir_array_temp_SSAValue2_28_1) / (SSAValue25pp2);
			parallel_ir_array_temp_SSAValue2_30_2 = SSAValue12;
			SSAValue20.ARRAYELEM(parfor_index_1_27) = parallel_ir_array_temp_SSAValue2_30_2;
		}
		;
		d = SSAValue20;
		parallel_ir_save_array_len_1_31 = d.ARRAYSIZE(1);
		parallel_ir_reduction_output_31 = DBL_MAX;
		for ( parfor_index_1_31 = 1; parfor_index_1_31 <= (int64_t)parallel_ir_save_array_len_1_31; parfor_index_1_31 += 1)
		{
			;
			parallel_ir_array_temp__4_33_1 = d.ARRAYELEM(parfor_index_1_31);
			parallel_ir_reduction_output_31 = std::min((double)parallel_ir_reduction_output_31,(double)parallel_ir_array_temp__4_33_1);
		}
		;
		m = parallel_ir_reduction_output_31;
		SSAValue23 = (m) - (SSAValue22);
		parallel_ir_save_array_len_1_35 = d.ARRAYSIZE(1);
		for ( parfor_index_1_35 = 1; parfor_index_1_35 <= (int64_t)parallel_ir_save_array_len_1_35; parfor_index_1_35 += 1)
		{
			;
			parallel_ir_array_temp__4_36_1 = d.ARRAYELEM(parfor_index_1_35);
			SSAValue15 = (parallel_ir_array_temp__4_36_1) - (m);
			parallel_ir_array_temp__4_38_2 = SSAValue15;
			d.ARRAYELEM(parfor_index_1_35) = parallel_ir_array_temp__4_38_2;
		}
		;
		SSAValue3 = d;
		parallel_ir_save_array_len_1_39 = SSAValue3.ARRAYSIZE(1);
		for ( parfor_index_1_39 = 1; parfor_index_1_39 <= (int64_t)parallel_ir_save_array_len_1_39; parfor_index_1_39 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue31_40_1 = SSAValue3.ARRAYELEM(parfor_index_1_39);
			SSAValue17 = exp(parallel_ir_array_temp_SSAValue31_40_1);
			parallel_ir_array_temp_SSAValue31_42_2 = SSAValue17;
			SSAValue3.ARRAYELEM(parfor_index_1_39) = parallel_ir_array_temp_SSAValue31_42_2;
		}
		;
		SSAValue4 = SSAValue3;
		parallel_ir_save_array_len_1_43 = SSAValue4.ARRAYSIZE(1);
		parallel_ir_reduction_output_43 = 0.0;
		for ( parfor_index_1_43 = 1; parfor_index_1_43 <= (int64_t)parallel_ir_save_array_len_1_43; parfor_index_1_43 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue32_45_1 = SSAValue4.ARRAYELEM(parfor_index_1_43);
			parallel_ir_reduction_output_43 = (parallel_ir_reduction_output_43+parallel_ir_array_temp_SSAValue32_45_1);
		}
		;
		SSAValue5 = parallel_ir_reduction_output_43;
		SSAValue24 = log(SSAValue5);;
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
	int64_t SSAValue9;
	int64_t SSAValue11;
	bool SSAValue13;
	int64_t SSAValue14;
	int64_t SSAValue16;
	int64_t SSAValue18;
	int64_t parfor_index_1_1;
	int64_t parallel_ir_save_array_len_1_1;
	double SSAValue19;
	double parallel_ir_array_temp__10_4_2;
	int64_t ppip273p280;
	double parallel_ir_reduction_input_5_1;
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
	int32_t SSAValue28;
	double SSAValue0;
	double SSAValue25pp2;
	double SSAValue1;
	double SSAValue2;
	j2c_array< double >  SSAValue3;
	j2c_array< double >  SSAValue4;
	double SSAValue5;
	double parallel_ir_array_temp__7_16_1;
	int64_t parfor_index_1_15;
	int64_t parallel_ir_save_array_len_1_15;
	double SSAValue6;
	j2c_array< double >  parallel_ir_new_array_name_15_1;
	double parallel_ir_array_temp__23_18_1;
	double parallel_ir_array_temp_SSAValue17_20_1;
	int64_t parfor_index_1_19;
	int64_t parallel_ir_save_array_len_1_19;
	int32_t SSAValue7;
	double SSAValue8;
	double parallel_ir_array_temp_SSAValue17_22_2;
	double parallel_ir_array_temp_SSAValue18_24_1;
	int64_t parfor_index_1_23;
	int64_t parallel_ir_save_array_len_1_23;
	double SSAValue10;
	double parallel_ir_array_temp_SSAValue18_26_2;
	double parallel_ir_array_temp_SSAValue2_28_1;
	int64_t parfor_index_1_27;
	int64_t parallel_ir_save_array_len_1_27;
	double SSAValue12;
	double parallel_ir_array_temp_SSAValue2_30_2;
	int64_t parfor_index_1_31;
	double parallel_ir_array_temp__4_33_1;
	int64_t parallel_ir_save_array_len_1_31;
	double parallel_ir_reduction_output_31;
	double parallel_ir_array_temp__4_36_1;
	int64_t parfor_index_1_35;
	int64_t parallel_ir_save_array_len_1_35;
	double SSAValue15;
	double parallel_ir_array_temp__4_38_2;
	double parallel_ir_array_temp_SSAValue31_40_1;
	int64_t parfor_index_1_39;
	int64_t parallel_ir_save_array_len_1_39;
	double SSAValue17;
	double parallel_ir_array_temp_SSAValue31_42_2;
	int64_t parfor_index_1_43;
	double parallel_ir_array_temp_SSAValue32_45_1;
	int64_t parallel_ir_save_array_len_1_43;
	double parallel_ir_reduction_output_43;
	int32_t __hpat_num_pes;
	int32_t __hpat_node_id;
	int64_t __hpat_dist_arr_start_1;
	int64_t __hpat_dist_arr_div_1;
	int64_t __hpat_dist_arr_count_1;
	int64_t __hpat_loop_start_1;
	int64_t __hpat_loop_end_1;
	int64_t __hpat_loop_div_1;
	int64_t __hpat_loop_start_5;
	int64_t __hpat_loop_end_5;
	int64_t __hpat_loop_div_5;
	double __hpat_reduce_2;
	std::random_device cgen_rand_device;
	std::uniform_real_distribution<double> cgen_distribution(0.0,1.0);
	std::normal_distribution<double> cgen_n_distribution(0.0,1.0);
	std::default_random_engine cgen_rand_generator(cgen_rand_device());
	;;
	MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
	MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
	SSAValue14 = (1) - (1);
	exps = 0.0;
	__hpat_dist_arr_div_1 = (n) / (__hpat_num_pes);
	__hpat_dist_arr_start_1 = (__hpat_node_id) * (__hpat_dist_arr_div_1);
	__hpat_dist_arr_count_1 = ((__hpat_node_id==__hpat_num_pes-1) ? n-__hpat_node_id*__hpat_dist_arr_div_1 : __hpat_dist_arr_div_1);
	arr = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
	parallel_ir_save_array_len_1_1 = n;
	__hpat_loop_div_1 = (parallel_ir_save_array_len_1_1) / (__hpat_num_pes);
	__hpat_loop_start_1 = ((__hpat_node_id) * (__hpat_loop_div_1)) + (1);
	__hpat_loop_end_1 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_1_1 : (__hpat_node_id+1)*__hpat_loop_div_1);
	for ( parfor_index_1_1 = __hpat_loop_start_1; parfor_index_1_1 <= (int64_t)__hpat_loop_end_1; parfor_index_1_1 += 1)
	{
		;
		SSAValue19 = cgen_distribution(cgen_rand_generator);
		;
		parallel_ir_array_temp__10_4_2 = SSAValue19;
		arr.ARRAYELEM(((parfor_index_1_1) - (__hpat_loop_start_1)) + (1)) = parallel_ir_array_temp__10_4_2;
	}
	;
	X = arr;
	points = _Base_vect(-1.0,2.0,5.0);
	N = points.ARRAYSIZE(1);
	b = 0.5;
	SSAValue9 = n;
	SSAValue13 = (1) <= (SSAValue9);
	SSAValue11 = (SSAValue13) ? (SSAValue9) : (SSAValue14);
	SSAValue16 = (SSAValue11) - (1);
	SSAValue18 = (SSAValue16) + (1);
	__hpat_loop_div_5 = (SSAValue18) / (__hpat_num_pes);
	__hpat_loop_start_5 = ((__hpat_node_id) * (__hpat_loop_div_5)) + (1);
	__hpat_loop_end_5 = ((__hpat_node_id==__hpat_num_pes-1) ? SSAValue18 : (__hpat_node_id+1)*__hpat_loop_div_5);
	#pragma simd reduction(+: exps)
	for ( ppip273p280 = __hpat_loop_start_5; ppip273p280 <= (int64_t)__hpat_loop_end_5; ppip273p280 += 1)
	{
		;
		SSAValue28 = (int32_t)(2);
		SSAValue21 = pow(b, SSAValue28);
		SSAValue0 = (double)2;
		SSAValue25pp2 = (SSAValue0) * (SSAValue21);
		SSAValue1 = (double)N;
		SSAValue2 = (b) * (SSAValue1);
		SSAValue22 = log(SSAValue2);;
		SSAValue16pp1 = X.ARRAYELEM(((ppip273p280) - (__hpat_loop_start_5)) + (1));
		parallel_ir_save_array_len_1_15 = points.ARRAYSIZE(1);
		parallel_ir_new_array_name_15_1 = j2c_array<double>::new_j2c_array_1d(NULL, parallel_ir_save_array_len_1_15);
		for ( parfor_index_1_15 = 1; parfor_index_1_15 <= (int64_t)parallel_ir_save_array_len_1_15; parfor_index_1_15 += 1)
		{
			;
			parallel_ir_array_temp__7_16_1 = points.ARRAYELEM(parfor_index_1_15);
			SSAValue6 = (SSAValue16pp1) - (parallel_ir_array_temp__7_16_1);
			parallel_ir_array_temp__23_18_1 = SSAValue6;
			parallel_ir_new_array_name_15_1.ARRAYELEM(parfor_index_1_15) = parallel_ir_array_temp__23_18_1;
		}
		;
		SSAValue26 = parallel_ir_new_array_name_15_1;
		parallel_ir_save_array_len_1_19 = SSAValue26.ARRAYSIZE(1);
		for ( parfor_index_1_19 = 1; parfor_index_1_19 <= (int64_t)parallel_ir_save_array_len_1_19; parfor_index_1_19 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue17_20_1 = SSAValue26.ARRAYELEM(parfor_index_1_19);
			SSAValue7 = (int32_t)(2);
			SSAValue8 = pow(parallel_ir_array_temp_SSAValue17_20_1, SSAValue7);
			parallel_ir_array_temp_SSAValue17_22_2 = SSAValue8;
			SSAValue26.ARRAYELEM(parfor_index_1_19) = parallel_ir_array_temp_SSAValue17_22_2;
		}
		;
		SSAValue27 = SSAValue26;
		parallel_ir_save_array_len_1_23 = SSAValue27.ARRAYSIZE(1);
		for ( parfor_index_1_23 = 1; parfor_index_1_23 <= (int64_t)parallel_ir_save_array_len_1_23; parfor_index_1_23 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue18_24_1 = SSAValue27.ARRAYELEM(parfor_index_1_23);
			SSAValue10 = -(parallel_ir_array_temp_SSAValue18_24_1);
			parallel_ir_array_temp_SSAValue18_26_2 = SSAValue10;
			SSAValue27.ARRAYELEM(parfor_index_1_23) = parallel_ir_array_temp_SSAValue18_26_2;
		}
		;
		SSAValue20 = SSAValue27;
		parallel_ir_save_array_len_1_27 = SSAValue20.ARRAYSIZE(1);
		for ( parfor_index_1_27 = 1; parfor_index_1_27 <= (int64_t)parallel_ir_save_array_len_1_27; parfor_index_1_27 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue2_28_1 = SSAValue20.ARRAYELEM(parfor_index_1_27);
			SSAValue12 = (parallel_ir_array_temp_SSAValue2_28_1) / (SSAValue25pp2);
			parallel_ir_array_temp_SSAValue2_30_2 = SSAValue12;
			SSAValue20.ARRAYELEM(parfor_index_1_27) = parallel_ir_array_temp_SSAValue2_30_2;
		}
		;
		d = SSAValue20;
		parallel_ir_save_array_len_1_31 = d.ARRAYSIZE(1);
		parallel_ir_reduction_output_31 = DBL_MAX;
		for ( parfor_index_1_31 = 1; parfor_index_1_31 <= (int64_t)parallel_ir_save_array_len_1_31; parfor_index_1_31 += 1)
		{
			;
			parallel_ir_array_temp__4_33_1 = d.ARRAYELEM(parfor_index_1_31);
			parallel_ir_reduction_output_31 = std::min((double)parallel_ir_reduction_output_31,(double)parallel_ir_array_temp__4_33_1);
		}
		;
		m = parallel_ir_reduction_output_31;
		SSAValue23 = (m) - (SSAValue22);
		parallel_ir_save_array_len_1_35 = d.ARRAYSIZE(1);
		for ( parfor_index_1_35 = 1; parfor_index_1_35 <= (int64_t)parallel_ir_save_array_len_1_35; parfor_index_1_35 += 1)
		{
			;
			parallel_ir_array_temp__4_36_1 = d.ARRAYELEM(parfor_index_1_35);
			SSAValue15 = (parallel_ir_array_temp__4_36_1) - (m);
			parallel_ir_array_temp__4_38_2 = SSAValue15;
			d.ARRAYELEM(parfor_index_1_35) = parallel_ir_array_temp__4_38_2;
		}
		;
		SSAValue3 = d;
		parallel_ir_save_array_len_1_39 = SSAValue3.ARRAYSIZE(1);
		for ( parfor_index_1_39 = 1; parfor_index_1_39 <= (int64_t)parallel_ir_save_array_len_1_39; parfor_index_1_39 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue31_40_1 = SSAValue3.ARRAYELEM(parfor_index_1_39);
			SSAValue17 = exp(parallel_ir_array_temp_SSAValue31_40_1);
			parallel_ir_array_temp_SSAValue31_42_2 = SSAValue17;
			SSAValue3.ARRAYELEM(parfor_index_1_39) = parallel_ir_array_temp_SSAValue31_42_2;
		}
		;
		SSAValue4 = SSAValue3;
		parallel_ir_save_array_len_1_43 = SSAValue4.ARRAYSIZE(1);
		parallel_ir_reduction_output_43 = 0.0;
		for ( parfor_index_1_43 = 1; parfor_index_1_43 <= (int64_t)parallel_ir_save_array_len_1_43; parfor_index_1_43 += 1)
		{
			;
			parallel_ir_array_temp_SSAValue32_45_1 = SSAValue4.ARRAYELEM(parfor_index_1_43);
			parallel_ir_reduction_output_43 = (parallel_ir_reduction_output_43+parallel_ir_array_temp_SSAValue32_45_1);
		}
		;
		SSAValue5 = parallel_ir_reduction_output_43;
		SSAValue24 = log(SSAValue5);;
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
