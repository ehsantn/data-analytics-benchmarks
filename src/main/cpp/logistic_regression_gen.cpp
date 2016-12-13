#include <random>
#include <mpi.h>
#include <mkl.h>
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
	j2c_array< double >  labels;
	j2c_array< double >  w;
	int64_t ptempp;
	j2c_array< double >  arr;
	j2c_array< double >  pparr_12p280;
	j2c_array< double >  pparr_14p282;
	j2c_array< double >  SSAValue10;
	j2c_array< double >  SSAValue15;
	j2c_array< double >  SSAValue21;
	int64_t SSAValue24;
	j2c_array< double >  SSAValue0;
	TupleInt64Int64 SSAValue1;
	TupleInt64Int64 SSAValue2;
	bool SSAValue3;
	int64_t SSAValue4;
	int64_t SSAValue5;
	bool SSAValue6;
	bool SSAValue7;
	j2c_array< double >  SSAValue8;
	int64_t SSAValue9;
	j2c_array< double >  SSAValue11;
	int64_t parfor_index_1_33;
	int64_t parallel_ir_save_array_len_1_33;
	double SSAValue12;
	double parallel_ir_array_temp__10_36_2;
	int64_t parfor_index_1_37;
	int64_t parfor_index_2_37;
	int64_t parallel_ir_save_array_len_1_37;
	int64_t parallel_ir_save_array_len_2_37;
	double SSAValue13;
	double parallel_ir_array_temp__12_40_2;
	int64_t parfor_index_1_41;
	int64_t parallel_ir_save_array_len_1_41;
	double SSAValue14;
	double parallel_ir_array_temp__14_44_2;
	double parallel_ir_array_temp__14_46_1;
	int64_t parallel_ir_save_array_len_1_45;
	double SSAValue16;
	double parallel_ir_array_temp__14_48_2;
	double parallel_ir_array_temp_SSAValue32_51_1;
	int64_t parallel_ir_save_array_len_1_50;
	double SSAValue17;
	double parallel_ir_array_temp_SSAValue32_53_2;
	double parallel_ir_array_temp__5_56_1;
	int64_t parfor_index_1_55;
	int64_t parfor_index_2_55;
	int64_t parallel_ir_save_array_len_1_55;
	int64_t parallel_ir_save_array_len_2_55;
	double SSAValue18;
	j2c_array< double >  parallel_ir_new_array_name_55_1;
	double parallel_ir_array_temp__83_58_1;
	double parallel_ir_array_temp_SSAValue16_60_1;
	double parallel_ir_array_temp_SSAValue45_61_1;
	double SSAValue19;
	double parallel_ir_array_temp_SSAValue16_63_2;
	double parallel_ir_array_temp_SSAValue46_66_1;
	double SSAValue20;
	double parallel_ir_array_temp_SSAValue46_68_2;
	double parallel_ir_array_temp_SSAValue47_71_1;
	double SSAValue22;
	double parallel_ir_array_temp_SSAValue47_73_2;
	double parallel_ir_array_temp_SSAValue48_76_1;
	double SSAValue23;
	double parallel_ir_array_temp_SSAValue48_78_2;
	double parallel_ir_array_temp_SSAValue49_81_1;
	double SSAValue25;
	double parallel_ir_array_temp_SSAValue49_83_2;
	double parallel_ir_array_temp_SSAValue50_86_1;
	double parallel_ir_array_temp__5_87_1;
	double SSAValue26;
	double parallel_ir_array_temp_SSAValue50_89_2;
	double parallel_ir_array_temp__7_92_1;
	double parallel_ir_array_temp_SSAValue53_93_1;
	int64_t parfor_index_1_91;
	int64_t parfor_index_2_91;
	int64_t parallel_ir_save_array_len_1_91;
	int64_t parallel_ir_save_array_len_2_91;
	double SSAValue27;
	j2c_array< double >  parallel_ir_new_array_name_91_1;
	double parallel_ir_array_temp__136_95_1;
	int32_t __hpat_num_pes;
	int32_t __hpat_node_id;
	int64_t __hpat_dist_arr_start_1;
	int64_t __hpat_dist_arr_div_1;
	int64_t __hpat_dist_arr_count_1;
	int64_t __hpat_loop_start_33;
	int64_t __hpat_loop_end_33;
	int64_t __hpat_loop_div_33;
	int64_t __hpat_dist_arr_start_2;
	int64_t __hpat_dist_arr_div_2;
	int64_t __hpat_dist_arr_count_2;
	TupleInt64Int64 __hpat_dist_tup_var_2;
	int64_t __hpat_dist_arr_start_3;
	int64_t __hpat_dist_arr_div_3;
	int64_t __hpat_dist_arr_count_3;
	int64_t __hpat_loop_start_37;
	int64_t __hpat_loop_end_37;
	int64_t __hpat_loop_div_37;
	int64_t __hpat_bcast_size_74;
	int64_t __hpat_dist_arr_start_4;
	int64_t __hpat_dist_arr_div_4;
	int64_t __hpat_dist_arr_count_4;
	int64_t __hpat_dist_arr_start_5;
	int64_t __hpat_dist_arr_div_5;
	int64_t __hpat_dist_arr_count_5;
	int64_t __hpat_loop_start_55;
	int64_t __hpat_loop_end_55;
	int64_t __hpat_loop_div_55;
	j2c_array< double >  __hpat_gemm_reduce_6;
	int64_t __hpat_gemm_reduce_size_6;
	std::random_device cgen_rand_device;
	std::uniform_real_distribution<double> cgen_distribution(0.0,1.0);
	std::normal_distribution<double> cgen_n_distribution(0.0,1.0);
	std::default_random_engine cgen_rand_generator(cgen_rand_device());
	;;
	MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
	MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
	SSAValue1 = TupleInt64Int64{1, N};
	SSAValue2 = TupleInt64Int64{1, 10};
	SSAValue3 = (1) <= (iterations);
	SSAValue4 = (1) - (1);
	SSAValue24 = (SSAValue3) ? (iterations) : (SSAValue4);
	SSAValue5 = (SSAValue24) + (1);
	__hpat_dist_arr_div_1 = (N) / (__hpat_num_pes);
	__hpat_dist_arr_start_1 = (__hpat_node_id) * (__hpat_dist_arr_div_1);
	__hpat_dist_arr_count_1 = ((__hpat_node_id==__hpat_num_pes-1) ? N-__hpat_node_id*__hpat_dist_arr_div_1 : __hpat_dist_arr_div_1);
	arr = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
	__hpat_loop_div_33 = (N) / (__hpat_num_pes);
	__hpat_loop_start_33 = ((__hpat_node_id) * (__hpat_loop_div_33)) + (1);
	__hpat_loop_end_33 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_33);
	for ( parfor_index_1_33 = __hpat_loop_start_33; parfor_index_1_33 <= (int64_t)__hpat_loop_end_33; parfor_index_1_33 += 1)
	{
		;
		SSAValue12 = cgen_distribution(cgen_rand_generator);
		;
		arr.ARRAYELEM(((parfor_index_1_33) - (__hpat_loop_start_33)) + (1)) = SSAValue12;
	}
	;
	__hpat_dist_arr_div_2 = (N) / (__hpat_num_pes);
	__hpat_dist_arr_start_2 = (__hpat_node_id) * (__hpat_dist_arr_div_2);
	__hpat_dist_arr_count_2 = ((__hpat_node_id==__hpat_num_pes-1) ? N-__hpat_node_id*__hpat_dist_arr_div_2 : __hpat_dist_arr_div_2);
	__hpat_dist_tup_var_2 = TupleInt64Int64{1, __hpat_dist_arr_count_2};
	labels = arr.reshape(__hpat_dist_tup_var_2.f0,__hpat_dist_tup_var_2.f1);
	;
	__hpat_dist_arr_div_3 = (N) / (__hpat_num_pes);
	__hpat_dist_arr_start_3 = (__hpat_node_id) * (__hpat_dist_arr_div_3);
	__hpat_dist_arr_count_3 = ((__hpat_node_id==__hpat_num_pes-1) ? N-__hpat_node_id*__hpat_dist_arr_div_3 : __hpat_dist_arr_div_3);
	pparr_12p280 = j2c_array<double>::new_j2c_array_2d(NULL, 10, __hpat_dist_arr_count_3);
	pparr_14p282 = j2c_array<double>::new_j2c_array_1d(NULL, 10);
	__hpat_loop_div_37 = (N) / (__hpat_num_pes);
	__hpat_loop_start_37 = ((__hpat_node_id) * (__hpat_loop_div_37)) + (1);
	__hpat_loop_end_37 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_37);
	for ( parfor_index_2_37 = __hpat_loop_start_37; parfor_index_2_37 <= (int64_t)__hpat_loop_end_37; parfor_index_2_37 += 1)
	{
		for ( parfor_index_1_37 = 1; parfor_index_1_37 <= (int64_t)10; parfor_index_1_37 += 1)
		{
			;
			SSAValue13 = cgen_distribution(cgen_rand_generator);
			;
			pparr_12p280.ARRAYELEM(parfor_index_1_37,((parfor_index_2_37) - (__hpat_loop_start_37)) + (1)) = SSAValue13;
		}
	}
	;
	if (!((__hpat_node_id == 0))) goto label74;
	for ( parfor_index_1_41 = 1; parfor_index_1_41 <= (int64_t)10; parfor_index_1_41 += 1)
	{
		;
		SSAValue14 = cgen_distribution(cgen_rand_generator);
		;
		SSAValue16 = (2.0) * (SSAValue14);
		SSAValue17 = (SSAValue16) - (1.0);
		pparr_14p282.ARRAYELEM(parfor_index_1_41) = SSAValue17;
	}
	;
	label74 : ;
	__hpat_bcast_size_74 = 10;
	MPI_Bcast(pparr_14p282.data, __hpat_bcast_size_74, MPI_DOUBLE, 0, MPI_COMM_WORLD);;
	w = pparr_14p282.reshape(SSAValue2.f0,SSAValue2.f1);
	;
	0; double __hpat_t1 = MPI_Wtime();
	ptempp = 1;
	label38 : ;
	SSAValue6 = (ptempp == SSAValue5);
	SSAValue7 = !(SSAValue6);
	if (!(SSAValue7)) goto label73;
	ptempp = (ptempp) + (1);
	__hpat_dist_arr_div_4 = (N) / (__hpat_num_pes);
	__hpat_dist_arr_start_4 = (__hpat_node_id) * (__hpat_dist_arr_div_4);
	__hpat_dist_arr_count_4 = ((__hpat_node_id==__hpat_num_pes-1) ? N-__hpat_node_id*__hpat_dist_arr_div_4 : __hpat_dist_arr_div_4);
	SSAValue15 = j2c_array<double>::new_j2c_array_2d(NULL, 1, __hpat_dist_arr_count_4);
	SSAValue8 = SSAValue15; cblas_dgemm((CBLAS_ORDER)102,(CBLAS_TRANSPOSE)111,(CBLAS_TRANSPOSE)111,w.ARRAYSIZE(1),pparr_12p280.ARRAYSIZE(2),w.ARRAYSIZE(2),1.0,
		w.data, w.ARRAYSIZE(1), pparr_12p280.data, pparr_12p280.ARRAYSIZE(1), 0.0, SSAValue15.data, w.ARRAYSIZE(1));
	parallel_ir_save_array_len_1_55 = 1;
	parallel_ir_save_array_len_2_55 = N;
	__hpat_dist_arr_div_5 = (parallel_ir_save_array_len_2_55) / (__hpat_num_pes);
	__hpat_dist_arr_start_5 = (__hpat_node_id) * (__hpat_dist_arr_div_5);
	__hpat_dist_arr_count_5 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_2_55-__hpat_node_id*__hpat_dist_arr_div_5 : __hpat_dist_arr_div_5);
	parallel_ir_new_array_name_55_1 = j2c_array<double>::new_j2c_array_2d(NULL, parallel_ir_save_array_len_1_55, __hpat_dist_arr_count_5);
	__hpat_loop_div_55 = (parallel_ir_save_array_len_2_55) / (__hpat_num_pes);
	__hpat_loop_start_55 = ((__hpat_node_id) * (__hpat_loop_div_55)) + (1);
	__hpat_loop_end_55 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_2_55 : (__hpat_node_id+1)*__hpat_loop_div_55);
	for ( parfor_index_2_55 = __hpat_loop_start_55; parfor_index_2_55 <= (int64_t)__hpat_loop_end_55; parfor_index_2_55 += 1)
	{
		for ( parfor_index_1_55 = 1; parfor_index_1_55 <= (int64_t)parallel_ir_save_array_len_1_55; parfor_index_1_55 += 1)
		{
			;
			parallel_ir_array_temp__5_56_1 = labels.ARRAYELEM(parfor_index_1_55,((parfor_index_2_55) - (__hpat_loop_start_55)) + (1));
			SSAValue18 = -(parallel_ir_array_temp__5_56_1);
			parallel_ir_array_temp_SSAValue45_61_1 = SSAValue8.ARRAYELEM(parfor_index_1_55,((parfor_index_2_55) - (__hpat_loop_start_55)) + (1));
			SSAValue19 = (SSAValue18) * (parallel_ir_array_temp_SSAValue45_61_1);
			SSAValue20 = exp(SSAValue19);
			SSAValue22 = (1.0) + (SSAValue20);
			SSAValue23 = (1.0) / (SSAValue22);
			SSAValue25 = (SSAValue23) - (1.0);
			parallel_ir_array_temp__5_87_1 = labels.ARRAYELEM(parfor_index_1_55,((parfor_index_2_55) - (__hpat_loop_start_55)) + (1));
			SSAValue26 = (SSAValue25) * (parallel_ir_array_temp__5_87_1);
			parallel_ir_new_array_name_55_1.ARRAYELEM(parfor_index_1_55,((parfor_index_2_55) - (__hpat_loop_start_55)) + (1)) = SSAValue26;
		}
	}
	;
	SSAValue9 = parallel_ir_save_array_len_1_55;
	SSAValue21 = j2c_array<double>::new_j2c_array_2d(NULL, SSAValue9, 10);
	__hpat_gemm_reduce_6 = j2c_array<double>::new_j2c_array_2d(NULL, SSAValue9, 10);
	SSAValue11 = __hpat_gemm_reduce_6; cblas_dgemm((CBLAS_ORDER)102,(CBLAS_TRANSPOSE)111,(CBLAS_TRANSPOSE)112,parallel_ir_new_array_name_55_1.ARRAYSIZE(1),pparr_12p280.ARRAYSIZE(1),parallel_ir_new_array_name_55_1.ARRAYSIZE(2),1.0,
		parallel_ir_new_array_name_55_1.data, parallel_ir_new_array_name_55_1.ARRAYSIZE(1), pparr_12p280.data, pparr_12p280.ARRAYSIZE(1), 0.0, __hpat_gemm_reduce_6.data, parallel_ir_new_array_name_55_1.ARRAYSIZE(1));
	__hpat_gemm_reduce_size_6 = (SSAValue9) * (10);
	MPI_Allreduce(__hpat_gemm_reduce_6.data, SSAValue21.data, __hpat_gemm_reduce_size_6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
	SSAValue11 = SSAValue21;
	parallel_ir_new_array_name_91_1 = j2c_array<double>::new_j2c_array_2d(NULL, 1, 10);
	for ( parfor_index_2_91 = 1; parfor_index_2_91 <= (int64_t)10; parfor_index_2_91 += 1)
	{
		for ( parfor_index_1_91 = 1; parfor_index_1_91 <= (int64_t)1; parfor_index_1_91 += 1)
		{
			;
			parallel_ir_array_temp__7_92_1 = w.ARRAYELEM(parfor_index_1_91,parfor_index_2_91);
			parallel_ir_array_temp_SSAValue53_93_1 = SSAValue11.ARRAYELEM(parfor_index_1_91,parfor_index_2_91);
			SSAValue27 = (parallel_ir_array_temp__7_92_1) - (parallel_ir_array_temp_SSAValue53_93_1);
			parallel_ir_new_array_name_91_1.ARRAYELEM(parfor_index_1_91,parfor_index_2_91) = SSAValue27;
		}
	}
	;
	w = parallel_ir_new_array_name_91_1;
	label71 : ;
	goto label38;
	label73 : ;
	0; if(__hpat_node_id==0) printf("exec time %lf\n", MPI_Wtime()-__hpat_t1);;
	*ret0 = w;
	return;

}


void pplogistic_regressionp271_unaliased(int64_t iterations, int64_t N,  j2c_array< double >  * __restrict ret0)
{
	ppplogistic_regressionp271 pselfp;
	j2c_array< double >  labels;
	j2c_array< double >  w;
	int64_t ptempp;
	j2c_array< double >  arr;
	j2c_array< double >  pparr_12p280;
	j2c_array< double >  pparr_14p282;
	j2c_array< double >  SSAValue10;
	j2c_array< double >  SSAValue15;
	j2c_array< double >  SSAValue21;
	int64_t SSAValue24;
	j2c_array< double >  SSAValue0;
	TupleInt64Int64 SSAValue1;
	TupleInt64Int64 SSAValue2;
	bool SSAValue3;
	int64_t SSAValue4;
	int64_t SSAValue5;
	bool SSAValue6;
	bool SSAValue7;
	j2c_array< double >  SSAValue8;
	int64_t SSAValue9;
	j2c_array< double >  SSAValue11;
	int64_t parfor_index_1_33;
	int64_t parallel_ir_save_array_len_1_33;
	double SSAValue12;
	double parallel_ir_array_temp__10_36_2;
	int64_t parfor_index_1_37;
	int64_t parfor_index_2_37;
	int64_t parallel_ir_save_array_len_1_37;
	int64_t parallel_ir_save_array_len_2_37;
	double SSAValue13;
	double parallel_ir_array_temp__12_40_2;
	int64_t parfor_index_1_41;
	int64_t parallel_ir_save_array_len_1_41;
	double SSAValue14;
	double parallel_ir_array_temp__14_44_2;
	double parallel_ir_array_temp__14_46_1;
	int64_t parallel_ir_save_array_len_1_45;
	double SSAValue16;
	double parallel_ir_array_temp__14_48_2;
	double parallel_ir_array_temp_SSAValue32_51_1;
	int64_t parallel_ir_save_array_len_1_50;
	double SSAValue17;
	double parallel_ir_array_temp_SSAValue32_53_2;
	double parallel_ir_array_temp__5_56_1;
	int64_t parfor_index_1_55;
	int64_t parfor_index_2_55;
	int64_t parallel_ir_save_array_len_1_55;
	int64_t parallel_ir_save_array_len_2_55;
	double SSAValue18;
	j2c_array< double >  parallel_ir_new_array_name_55_1;
	double parallel_ir_array_temp__83_58_1;
	double parallel_ir_array_temp_SSAValue16_60_1;
	double parallel_ir_array_temp_SSAValue45_61_1;
	double SSAValue19;
	double parallel_ir_array_temp_SSAValue16_63_2;
	double parallel_ir_array_temp_SSAValue46_66_1;
	double SSAValue20;
	double parallel_ir_array_temp_SSAValue46_68_2;
	double parallel_ir_array_temp_SSAValue47_71_1;
	double SSAValue22;
	double parallel_ir_array_temp_SSAValue47_73_2;
	double parallel_ir_array_temp_SSAValue48_76_1;
	double SSAValue23;
	double parallel_ir_array_temp_SSAValue48_78_2;
	double parallel_ir_array_temp_SSAValue49_81_1;
	double SSAValue25;
	double parallel_ir_array_temp_SSAValue49_83_2;
	double parallel_ir_array_temp_SSAValue50_86_1;
	double parallel_ir_array_temp__5_87_1;
	double SSAValue26;
	double parallel_ir_array_temp_SSAValue50_89_2;
	double parallel_ir_array_temp__7_92_1;
	double parallel_ir_array_temp_SSAValue53_93_1;
	int64_t parfor_index_1_91;
	int64_t parfor_index_2_91;
	int64_t parallel_ir_save_array_len_1_91;
	int64_t parallel_ir_save_array_len_2_91;
	double SSAValue27;
	j2c_array< double >  parallel_ir_new_array_name_91_1;
	double parallel_ir_array_temp__136_95_1;
	int32_t __hpat_num_pes;
	int32_t __hpat_node_id;
	int64_t __hpat_dist_arr_start_1;
	int64_t __hpat_dist_arr_div_1;
	int64_t __hpat_dist_arr_count_1;
	int64_t __hpat_loop_start_33;
	int64_t __hpat_loop_end_33;
	int64_t __hpat_loop_div_33;
	int64_t __hpat_dist_arr_start_2;
	int64_t __hpat_dist_arr_div_2;
	int64_t __hpat_dist_arr_count_2;
	TupleInt64Int64 __hpat_dist_tup_var_2;
	int64_t __hpat_dist_arr_start_3;
	int64_t __hpat_dist_arr_div_3;
	int64_t __hpat_dist_arr_count_3;
	int64_t __hpat_loop_start_37;
	int64_t __hpat_loop_end_37;
	int64_t __hpat_loop_div_37;
	int64_t __hpat_bcast_size_74;
	int64_t __hpat_dist_arr_start_4;
	int64_t __hpat_dist_arr_div_4;
	int64_t __hpat_dist_arr_count_4;
	int64_t __hpat_dist_arr_start_5;
	int64_t __hpat_dist_arr_div_5;
	int64_t __hpat_dist_arr_count_5;
	int64_t __hpat_loop_start_55;
	int64_t __hpat_loop_end_55;
	int64_t __hpat_loop_div_55;
	j2c_array< double >  __hpat_gemm_reduce_6;
	int64_t __hpat_gemm_reduce_size_6;
	std::random_device cgen_rand_device;
	std::uniform_real_distribution<double> cgen_distribution(0.0,1.0);
	std::normal_distribution<double> cgen_n_distribution(0.0,1.0);
	std::default_random_engine cgen_rand_generator(cgen_rand_device());
	;;
	MPI_Comm_size(MPI_COMM_WORLD,&__hpat_num_pes);;
	MPI_Comm_rank(MPI_COMM_WORLD,&__hpat_node_id);;
	SSAValue1 = TupleInt64Int64{1, N};
	SSAValue2 = TupleInt64Int64{1, 10};
	SSAValue3 = (1) <= (iterations);
	SSAValue4 = (1) - (1);
	SSAValue24 = (SSAValue3) ? (iterations) : (SSAValue4);
	SSAValue5 = (SSAValue24) + (1);
	__hpat_dist_arr_div_1 = (N) / (__hpat_num_pes);
	__hpat_dist_arr_start_1 = (__hpat_node_id) * (__hpat_dist_arr_div_1);
	__hpat_dist_arr_count_1 = ((__hpat_node_id==__hpat_num_pes-1) ? N-__hpat_node_id*__hpat_dist_arr_div_1 : __hpat_dist_arr_div_1);
	arr = j2c_array<double>::new_j2c_array_1d(NULL, __hpat_dist_arr_count_1);
	__hpat_loop_div_33 = (N) / (__hpat_num_pes);
	__hpat_loop_start_33 = ((__hpat_node_id) * (__hpat_loop_div_33)) + (1);
	__hpat_loop_end_33 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_33);
	for ( parfor_index_1_33 = __hpat_loop_start_33; parfor_index_1_33 <= (int64_t)__hpat_loop_end_33; parfor_index_1_33 += 1)
	{
		;
		SSAValue12 = cgen_distribution(cgen_rand_generator);
		;
		arr.ARRAYELEM(((parfor_index_1_33) - (__hpat_loop_start_33)) + (1)) = SSAValue12;
	}
	;
	__hpat_dist_arr_div_2 = (N) / (__hpat_num_pes);
	__hpat_dist_arr_start_2 = (__hpat_node_id) * (__hpat_dist_arr_div_2);
	__hpat_dist_arr_count_2 = ((__hpat_node_id==__hpat_num_pes-1) ? N-__hpat_node_id*__hpat_dist_arr_div_2 : __hpat_dist_arr_div_2);
	__hpat_dist_tup_var_2 = TupleInt64Int64{1, __hpat_dist_arr_count_2};
	labels = arr.reshape(__hpat_dist_tup_var_2.f0,__hpat_dist_tup_var_2.f1);
	;
	__hpat_dist_arr_div_3 = (N) / (__hpat_num_pes);
	__hpat_dist_arr_start_3 = (__hpat_node_id) * (__hpat_dist_arr_div_3);
	__hpat_dist_arr_count_3 = ((__hpat_node_id==__hpat_num_pes-1) ? N-__hpat_node_id*__hpat_dist_arr_div_3 : __hpat_dist_arr_div_3);
	pparr_12p280 = j2c_array<double>::new_j2c_array_2d(NULL, 10, __hpat_dist_arr_count_3);
	pparr_14p282 = j2c_array<double>::new_j2c_array_1d(NULL, 10);
	__hpat_loop_div_37 = (N) / (__hpat_num_pes);
	__hpat_loop_start_37 = ((__hpat_node_id) * (__hpat_loop_div_37)) + (1);
	__hpat_loop_end_37 = ((__hpat_node_id==__hpat_num_pes-1) ? N : (__hpat_node_id+1)*__hpat_loop_div_37);
	for ( parfor_index_2_37 = __hpat_loop_start_37; parfor_index_2_37 <= (int64_t)__hpat_loop_end_37; parfor_index_2_37 += 1)
	{
		for ( parfor_index_1_37 = 1; parfor_index_1_37 <= (int64_t)10; parfor_index_1_37 += 1)
		{
			;
			SSAValue13 = cgen_distribution(cgen_rand_generator);
			;
			pparr_12p280.ARRAYELEM(parfor_index_1_37,((parfor_index_2_37) - (__hpat_loop_start_37)) + (1)) = SSAValue13;
		}
	}
	;
	if (!((__hpat_node_id == 0))) goto label74;
	for ( parfor_index_1_41 = 1; parfor_index_1_41 <= (int64_t)10; parfor_index_1_41 += 1)
	{
		;
		SSAValue14 = cgen_distribution(cgen_rand_generator);
		;
		SSAValue16 = (2.0) * (SSAValue14);
		SSAValue17 = (SSAValue16) - (1.0);
		pparr_14p282.ARRAYELEM(parfor_index_1_41) = SSAValue17;
	}
	;
	label74 : ;
	__hpat_bcast_size_74 = 10;
	MPI_Bcast(pparr_14p282.data, __hpat_bcast_size_74, MPI_DOUBLE, 0, MPI_COMM_WORLD);;
	w = pparr_14p282.reshape(SSAValue2.f0,SSAValue2.f1);
	;
	0; double __hpat_t1 = MPI_Wtime();
	ptempp = 1;
	label38 : ;
	SSAValue6 = (ptempp == SSAValue5);
	SSAValue7 = !(SSAValue6);
	if (!(SSAValue7)) goto label73;
	ptempp = (ptempp) + (1);
	__hpat_dist_arr_div_4 = (N) / (__hpat_num_pes);
	__hpat_dist_arr_start_4 = (__hpat_node_id) * (__hpat_dist_arr_div_4);
	__hpat_dist_arr_count_4 = ((__hpat_node_id==__hpat_num_pes-1) ? N-__hpat_node_id*__hpat_dist_arr_div_4 : __hpat_dist_arr_div_4);
	SSAValue15 = j2c_array<double>::new_j2c_array_2d(NULL, 1, __hpat_dist_arr_count_4);
	SSAValue8 = SSAValue15; cblas_dgemm((CBLAS_ORDER)102,(CBLAS_TRANSPOSE)111,(CBLAS_TRANSPOSE)111,w.ARRAYSIZE(1),pparr_12p280.ARRAYSIZE(2),w.ARRAYSIZE(2),1.0,
		w.data, w.ARRAYSIZE(1), pparr_12p280.data, pparr_12p280.ARRAYSIZE(1), 0.0, SSAValue15.data, w.ARRAYSIZE(1));
	parallel_ir_save_array_len_1_55 = 1;
	parallel_ir_save_array_len_2_55 = N;
	__hpat_dist_arr_div_5 = (parallel_ir_save_array_len_2_55) / (__hpat_num_pes);
	__hpat_dist_arr_start_5 = (__hpat_node_id) * (__hpat_dist_arr_div_5);
	__hpat_dist_arr_count_5 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_2_55-__hpat_node_id*__hpat_dist_arr_div_5 : __hpat_dist_arr_div_5);
	parallel_ir_new_array_name_55_1 = j2c_array<double>::new_j2c_array_2d(NULL, parallel_ir_save_array_len_1_55, __hpat_dist_arr_count_5);
	__hpat_loop_div_55 = (parallel_ir_save_array_len_2_55) / (__hpat_num_pes);
	__hpat_loop_start_55 = ((__hpat_node_id) * (__hpat_loop_div_55)) + (1);
	__hpat_loop_end_55 = ((__hpat_node_id==__hpat_num_pes-1) ? parallel_ir_save_array_len_2_55 : (__hpat_node_id+1)*__hpat_loop_div_55);
	for ( parfor_index_2_55 = __hpat_loop_start_55; parfor_index_2_55 <= (int64_t)__hpat_loop_end_55; parfor_index_2_55 += 1)
	{
		for ( parfor_index_1_55 = 1; parfor_index_1_55 <= (int64_t)parallel_ir_save_array_len_1_55; parfor_index_1_55 += 1)
		{
			;
			parallel_ir_array_temp__5_56_1 = labels.ARRAYELEM(parfor_index_1_55,((parfor_index_2_55) - (__hpat_loop_start_55)) + (1));
			SSAValue18 = -(parallel_ir_array_temp__5_56_1);
			parallel_ir_array_temp_SSAValue45_61_1 = SSAValue8.ARRAYELEM(parfor_index_1_55,((parfor_index_2_55) - (__hpat_loop_start_55)) + (1));
			SSAValue19 = (SSAValue18) * (parallel_ir_array_temp_SSAValue45_61_1);
			SSAValue20 = exp(SSAValue19);
			SSAValue22 = (1.0) + (SSAValue20);
			SSAValue23 = (1.0) / (SSAValue22);
			SSAValue25 = (SSAValue23) - (1.0);
			parallel_ir_array_temp__5_87_1 = labels.ARRAYELEM(parfor_index_1_55,((parfor_index_2_55) - (__hpat_loop_start_55)) + (1));
			SSAValue26 = (SSAValue25) * (parallel_ir_array_temp__5_87_1);
			parallel_ir_new_array_name_55_1.ARRAYELEM(parfor_index_1_55,((parfor_index_2_55) - (__hpat_loop_start_55)) + (1)) = SSAValue26;
		}
	}
	;
	SSAValue9 = parallel_ir_save_array_len_1_55;
	SSAValue21 = j2c_array<double>::new_j2c_array_2d(NULL, SSAValue9, 10);
	__hpat_gemm_reduce_6 = j2c_array<double>::new_j2c_array_2d(NULL, SSAValue9, 10);
	SSAValue11 = __hpat_gemm_reduce_6; cblas_dgemm((CBLAS_ORDER)102,(CBLAS_TRANSPOSE)111,(CBLAS_TRANSPOSE)112,parallel_ir_new_array_name_55_1.ARRAYSIZE(1),pparr_12p280.ARRAYSIZE(1),parallel_ir_new_array_name_55_1.ARRAYSIZE(2),1.0,
		parallel_ir_new_array_name_55_1.data, parallel_ir_new_array_name_55_1.ARRAYSIZE(1), pparr_12p280.data, pparr_12p280.ARRAYSIZE(1), 0.0, __hpat_gemm_reduce_6.data, parallel_ir_new_array_name_55_1.ARRAYSIZE(1));
	__hpat_gemm_reduce_size_6 = (SSAValue9) * (10);
	MPI_Allreduce(__hpat_gemm_reduce_6.data, SSAValue21.data, __hpat_gemm_reduce_size_6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);;
	SSAValue11 = SSAValue21;
	parallel_ir_new_array_name_91_1 = j2c_array<double>::new_j2c_array_2d(NULL, 1, 10);
	for ( parfor_index_2_91 = 1; parfor_index_2_91 <= (int64_t)10; parfor_index_2_91 += 1)
	{
		for ( parfor_index_1_91 = 1; parfor_index_1_91 <= (int64_t)1; parfor_index_1_91 += 1)
		{
			;
			parallel_ir_array_temp__7_92_1 = w.ARRAYELEM(parfor_index_1_91,parfor_index_2_91);
			parallel_ir_array_temp_SSAValue53_93_1 = SSAValue11.ARRAYELEM(parfor_index_1_91,parfor_index_2_91);
			SSAValue27 = (parallel_ir_array_temp__7_92_1) - (parallel_ir_array_temp_SSAValue53_93_1);
			parallel_ir_new_array_name_91_1.ARRAYELEM(parfor_index_1_91,parfor_index_2_91) = SSAValue27;
		}
	}
	;
	w = parallel_ir_new_array_name_91_1;
	label71 : ;
	goto label38;
	label73 : ;
	0; if(__hpat_node_id==0) printf("exec time %lf\n", MPI_Wtime()-__hpat_t1);;
	*ret0 = w;
	return;

}


extern "C" void _pplogistic_regressionp271_unaliased_(int run_where, int64_t iterations, int64_t N ,  j2c_array< double > ** __restrict ret0 , bool genMain = true)
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
		mainFile << "    _pplogistic_regressionp271_unaliased_(runwhere, iterations, N, &ret0, false);" << std::endl;
		mainFile << "    MPI_Finalize();" << std::endl;
		mainFile << "    return 0;" << std::endl;
		mainFile << "}" << std::endl;
		mainFile.close();
		std::ofstream mainFileSh(newMainSh.str());
		mainFileSh << "#!/bin/sh" << std::endl;
		mainFileSh << "mpiicpc -O3 -std=c++11 -I/usr/local/hdf5/include  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << " -mkl  -lm " << std::endl;
		mainFileSh.close();
	}

	*ret0 = new  j2c_array< double > ();

	pplogistic_regressionp271_unaliased(iterations, N, *ret0);
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
		mainFileSh << "mpiicpc -O3 -std=c++11 -I/usr/local/hdf5/include  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << " -mkl  -lm " << std::endl;
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
