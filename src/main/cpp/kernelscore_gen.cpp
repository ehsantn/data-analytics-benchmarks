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
    int64_t ppptempp_4p297;
    int64_t ppptempp_5p298;
    int64_t i;
    int64_t ppptempp_7p299;
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
    ppptempp_4p297 = 1;
    ppptempp_5p298 = 0;
    label8 : ;
    if (!(!((ppptempp_5p298 == SSAValue1)))) goto label28;
    SSAValue3 = (ppptempp_5p298) + (1);
    ppptempp_5p298 = SSAValue3;
    SSAValue11 = ppptempp_4p297;
    SSAValue12 = (ppptempp_4p297) + (1);
    ppptempp_7p299 = 1;
    SSAValue13 = (1) + (1);
    i = SSAValue11;
    ppptempp_7p299 = SSAValue13;
    SSAValue14 = (2) + (1);
    ppptempp_4p297 = SSAValue12;
    ppptempp_7p299 = SSAValue14;
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
    int64_t parfor_index_1_3;
    int64_t parallel_ir_save_array_len_1_3;
    double SSAValue19;
    double parallel_ir_array_temp__10_6_2;
    int64_t ppip273p280;
    double parallel_ir_reduction_input_7_1;
    double pptemp_neutral_valp281;
    j2c_array< double >  d;
    j2c_array< double >  SSAValue20;
    double SSAValue21;
    double SSAValue22;
    double SSAValue23;
    double SSAValue24;
    double SSAValue25;
    double SSAValue26;
    double SSAValue16pp1;
    int32_t SSAValue27;
    double SSAValue28;
    double SSAValue25pp2;
    double SSAValue29;
    double SSAValue30;
    double SSAValue31pp3;
    double SSAValue31;
    double parallel_ir_array_temp__6_18_1;
    int64_t parfor_index_1_17;
    int64_t parallel_ir_save_array_len_1_17;
    double SSAValue32;
    j2c_array< double >  parallel_ir_new_array_name_17_1;
    double parallel_ir_array_temp__23_20_1;
    double parallel_ir_array_temp_SSAValue17_22_1;
    int32_t SSAValue33;
    double SSAValue34;
    double parallel_ir_array_temp_SSAValue17_24_2;
    double parallel_ir_array_temp_SSAValue18_27_1;
    double SSAValue35;
    double parallel_ir_array_temp_SSAValue18_29_2;
    double parallel_ir_array_temp_SSAValue2_32_1;
    int64_t parfor_index_1_31;
    int64_t parallel_ir_save_array_len_1_31;
    double SSAValue36;
    double parallel_ir_array_temp_SSAValue2_34_2;
    double parallel_ir_array_temp__4_37_1;
    double parallel_ir_reduction_output_35;
    int64_t SSAValue37;
    int64_t SSAValue38;
    bool SSAValue39;
    bool SSAValue40;
    bool SSAValue41;
    bool SSAValue42;
    bool SSAValue43;
    bool SSAValue44;
    bool SSAValue45;
    bool SSAValue46;
    double SSAValue47;
    double SSAValue48;
    double SSAValue49;
    int64_t parfor_index_1_40;
    double parallel_ir_array_temp__4_42_1;
    int64_t parallel_ir_save_array_len_1_40;
    double parallel_ir_reduction_output_40;
    int64_t SSAValue50;
    int64_t SSAValue51;
    bool SSAValue52;
    bool SSAValue0;
    bool SSAValue1;
    bool SSAValue2;
    bool SSAValue3;
    bool SSAValue4;
    bool SSAValue5;
    bool SSAValue6;
    double SSAValue7;
    double SSAValue8;
    double SSAValue10;
    double parallel_ir_array_temp__4_45_1;
    int64_t parfor_index_1_44;
    int64_t parallel_ir_save_array_len_1_44;
    double SSAValue12;
    double parallel_ir_array_temp__4_47_2;
    double parallel_ir_array_temp_SSAValue32_49_1;
    double SSAValue15;
    double parallel_ir_array_temp_SSAValue32_51_2;
    double parallel_ir_array_temp_SSAValue33_55_1;
    double parallel_ir_reduction_output_53;
    double SSAValue17;
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
    X = arr;
    t1 = MPI_Wtime();
    points = _Base_vect(-1.0,2.0,5.0);
    N = points.ARRAYSIZE(1);
    b = 0.5;
    exps = 0.0;
    SSAValue9 = n;
    SSAValue13 = (1) <= (SSAValue9);
    SSAValue14 = (1) - (1);
    SSAValue11 = (SSAValue13) ? (SSAValue9) : (SSAValue14);
    SSAValue16 = (SSAValue11) - (1);
    SSAValue18 = (SSAValue16) + (1);
    __hpat_loop_div_7 = (SSAValue18) / (__hpat_num_pes);
    __hpat_loop_start_7 = ((__hpat_node_id) * (__hpat_loop_div_7)) + (1);
    __hpat_loop_end_7 = ((__hpat_node_id==__hpat_num_pes-1) ? SSAValue18 : (__hpat_node_id+1)*__hpat_loop_div_7);
    for ( ppip273p280 = __hpat_loop_start_7; ppip273p280 <= (int64_t)__hpat_loop_end_7; ppip273p280 += 1)
    {
        ;
        SSAValue16pp1 = X.ARRAYELEM(((ppip273p280) - (__hpat_loop_start_7)) + (1));
        parallel_ir_save_array_len_1_17 = points.ARRAYSIZE(1);
        parallel_ir_new_array_name_17_1 = j2c_array<double>::new_j2c_array_1d(NULL, parallel_ir_save_array_len_1_17);
        for ( parfor_index_1_17 = 1; parfor_index_1_17 <= (int64_t)parallel_ir_save_array_len_1_17; parfor_index_1_17 += 1)
        {
            ;
            parallel_ir_array_temp__6_18_1 = points.ARRAYELEM(parfor_index_1_17);
            SSAValue32 = (SSAValue16pp1) - (parallel_ir_array_temp__6_18_1);
            parallel_ir_array_temp__23_20_1 = SSAValue32;
            parallel_ir_array_temp_SSAValue17_22_1 = parallel_ir_array_temp__23_20_1;
            SSAValue33 = (int32_t)(2);
            pow(parallel_ir_array_temp_SSAValue17_22_1, SSAValue33);
            SSAValue34 = pow(parallel_ir_array_temp_SSAValue17_22_1, SSAValue33);
            parallel_ir_array_temp_SSAValue17_24_2 = SSAValue34;
            parallel_ir_array_temp_SSAValue18_27_1 = parallel_ir_array_temp_SSAValue17_24_2;
            SSAValue35 = -(parallel_ir_array_temp_SSAValue18_27_1);
            parallel_ir_array_temp_SSAValue18_29_2 = SSAValue35;
            parallel_ir_new_array_name_17_1.ARRAYELEM(parfor_index_1_17) = parallel_ir_array_temp_SSAValue18_29_2;
        }
        ;
        SSAValue20 = parallel_ir_new_array_name_17_1;
        SSAValue27 = (int32_t)(2);
        pow(b, SSAValue27);
        SSAValue22 = pow(b, SSAValue27);
        SSAValue28 = (double)2;
        SSAValue25pp2 = (SSAValue28) * (SSAValue22);
        parallel_ir_save_array_len_1_31 = SSAValue20.ARRAYSIZE(1);
        parallel_ir_reduction_output_35 = DBL_MAX;
        for ( parfor_index_1_31 = 1; parfor_index_1_31 <= (int64_t)parallel_ir_save_array_len_1_31; parfor_index_1_31 += 1)
        {
            ;
            parallel_ir_array_temp_SSAValue2_32_1 = SSAValue20.ARRAYELEM(parfor_index_1_31);
            SSAValue36 = (parallel_ir_array_temp_SSAValue2_32_1) / (SSAValue25pp2);
            parallel_ir_array_temp_SSAValue2_34_2 = SSAValue36;
            SSAValue20.ARRAYELEM(parfor_index_1_31) = parallel_ir_array_temp_SSAValue2_34_2;
            parallel_ir_array_temp__4_37_1 = parallel_ir_array_temp_SSAValue2_34_2;
            SSAValue37 = parallel_ir_array_temp__4_37_1;
            SSAValue38 = parallel_ir_reduction_output_35;
            SSAValue39 = (SSAValue38) < (0);
            SSAValue40 = (SSAValue37) < (0);
            SSAValue41 = !(SSAValue39);
            (SSAValue40) & (SSAValue41);
            SSAValue42 = (parallel_ir_array_temp__4_37_1) < (parallel_ir_reduction_output_35);
            SSAValue43 = (SSAValue40) & (SSAValue41);
            (SSAValue42) | (SSAValue43);
            SSAValue44 = (parallel_ir_array_temp__4_37_1) != (parallel_ir_array_temp__4_37_1);
            SSAValue45 = (parallel_ir_reduction_output_35) != (parallel_ir_reduction_output_35);
            SSAValue46 = (SSAValue42) | (SSAValue43);
            SSAValue47 = (SSAValue44) ? (parallel_ir_reduction_output_35) : (parallel_ir_array_temp__4_37_1);
            SSAValue48 = (SSAValue45) ? (parallel_ir_array_temp__4_37_1) : (parallel_ir_reduction_output_35);
            SSAValue49 = (SSAValue46) ? (SSAValue47) : (SSAValue48);
            parallel_ir_reduction_output_35 = SSAValue49;
        }
        ;
        SSAValue21 = parallel_ir_reduction_output_35;
        d = SSAValue20;
        SSAValue29 = (double)N;
        SSAValue30 = (b) * (SSAValue29);
        SSAValue23 = log(SSAValue30);;
        SSAValue24 = (SSAValue21) - (SSAValue23);
        parallel_ir_save_array_len_1_40 = d.ARRAYSIZE(1);
        parallel_ir_reduction_output_40 = DBL_MAX;
        for ( parfor_index_1_40 = 1; parfor_index_1_40 <= (int64_t)parallel_ir_save_array_len_1_40; parfor_index_1_40 += 1)
        {
            ;
            parallel_ir_array_temp__4_42_1 = d.ARRAYELEM(parfor_index_1_40);
            SSAValue50 = parallel_ir_array_temp__4_42_1;
            SSAValue51 = parallel_ir_reduction_output_40;
            SSAValue52 = (SSAValue51) < (0);
            SSAValue0 = (SSAValue50) < (0);
            SSAValue1 = !(SSAValue52);
            (SSAValue0) & (SSAValue1);
            SSAValue2 = (parallel_ir_array_temp__4_42_1) < (parallel_ir_reduction_output_40);
            SSAValue3 = (SSAValue0) & (SSAValue1);
            (SSAValue2) | (SSAValue3);
            SSAValue4 = (parallel_ir_array_temp__4_42_1) != (parallel_ir_array_temp__4_42_1);
            SSAValue5 = (parallel_ir_reduction_output_40) != (parallel_ir_reduction_output_40);
            SSAValue6 = (SSAValue2) | (SSAValue3);
            SSAValue7 = (SSAValue4) ? (parallel_ir_reduction_output_40) : (parallel_ir_array_temp__4_42_1);
            SSAValue8 = (SSAValue5) ? (parallel_ir_array_temp__4_42_1) : (parallel_ir_reduction_output_40);
            SSAValue10 = (SSAValue6) ? (SSAValue7) : (SSAValue8);
            parallel_ir_reduction_output_40 = SSAValue10;
        }
        ;
        SSAValue31pp3 = parallel_ir_reduction_output_40;
        parallel_ir_save_array_len_1_44 = d.ARRAYSIZE(1);
        parallel_ir_reduction_output_53 = 0.0;
        for ( parfor_index_1_44 = 1; parfor_index_1_44 <= (int64_t)parallel_ir_save_array_len_1_44; parfor_index_1_44 += 1)
        {
            ;
            parallel_ir_array_temp__4_45_1 = d.ARRAYELEM(parfor_index_1_44);
            SSAValue12 = (parallel_ir_array_temp__4_45_1) - (SSAValue31pp3);
            parallel_ir_array_temp__4_47_2 = SSAValue12;
            parallel_ir_array_temp_SSAValue32_49_1 = parallel_ir_array_temp__4_47_2;
            SSAValue15 = exp(parallel_ir_array_temp_SSAValue32_49_1);
            parallel_ir_array_temp_SSAValue32_51_2 = SSAValue15;
            parallel_ir_array_temp_SSAValue33_55_1 = parallel_ir_array_temp_SSAValue32_51_2;
            SSAValue17 = (parallel_ir_reduction_output_53) + (parallel_ir_array_temp_SSAValue33_55_1);
            parallel_ir_reduction_output_53 = SSAValue17;
        }
        ;
        SSAValue31 = parallel_ir_reduction_output_53;
        SSAValue25 = log(SSAValue31);;
        SSAValue26 = (SSAValue24) + (SSAValue25);
        exps = (exps) + (SSAValue26);
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
    int64_t parfor_index_1_3;
    int64_t parallel_ir_save_array_len_1_3;
    double SSAValue19;
    double parallel_ir_array_temp__10_6_2;
    int64_t ppip273p280;
    double parallel_ir_reduction_input_7_1;
    double pptemp_neutral_valp281;
    j2c_array< double >  d;
    j2c_array< double >  SSAValue20;
    double SSAValue21;
    double SSAValue22;
    double SSAValue23;
    double SSAValue24;
    double SSAValue25;
    double SSAValue26;
    double SSAValue16pp1;
    int32_t SSAValue27;
    double SSAValue28;
    double SSAValue25pp2;
    double SSAValue29;
    double SSAValue30;
    double SSAValue31pp3;
    double SSAValue31;
    double parallel_ir_array_temp__6_18_1;
    int64_t parfor_index_1_17;
    int64_t parallel_ir_save_array_len_1_17;
    double SSAValue32;
    j2c_array< double >  parallel_ir_new_array_name_17_1;
    double parallel_ir_array_temp__23_20_1;
    double parallel_ir_array_temp_SSAValue17_22_1;
    int32_t SSAValue33;
    double SSAValue34;
    double parallel_ir_array_temp_SSAValue17_24_2;
    double parallel_ir_array_temp_SSAValue18_27_1;
    double SSAValue35;
    double parallel_ir_array_temp_SSAValue18_29_2;
    double parallel_ir_array_temp_SSAValue2_32_1;
    int64_t parfor_index_1_31;
    int64_t parallel_ir_save_array_len_1_31;
    double SSAValue36;
    double parallel_ir_array_temp_SSAValue2_34_2;
    double parallel_ir_array_temp__4_37_1;
    double parallel_ir_reduction_output_35;
    int64_t SSAValue37;
    int64_t SSAValue38;
    bool SSAValue39;
    bool SSAValue40;
    bool SSAValue41;
    bool SSAValue42;
    bool SSAValue43;
    bool SSAValue44;
    bool SSAValue45;
    bool SSAValue46;
    double SSAValue47;
    double SSAValue48;
    double SSAValue49;
    int64_t parfor_index_1_40;
    double parallel_ir_array_temp__4_42_1;
    int64_t parallel_ir_save_array_len_1_40;
    double parallel_ir_reduction_output_40;
    int64_t SSAValue50;
    int64_t SSAValue51;
    bool SSAValue52;
    bool SSAValue0;
    bool SSAValue1;
    bool SSAValue2;
    bool SSAValue3;
    bool SSAValue4;
    bool SSAValue5;
    bool SSAValue6;
    double SSAValue7;
    double SSAValue8;
    double SSAValue10;
    double parallel_ir_array_temp__4_45_1;
    int64_t parfor_index_1_44;
    int64_t parallel_ir_save_array_len_1_44;
    double SSAValue12;
    double parallel_ir_array_temp__4_47_2;
    double parallel_ir_array_temp_SSAValue32_49_1;
    double SSAValue15;
    double parallel_ir_array_temp_SSAValue32_51_2;
    double parallel_ir_array_temp_SSAValue33_55_1;
    double parallel_ir_reduction_output_53;
    double SSAValue17;
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
    X = arr;
    points = _Base_vect(-1.0,2.0,5.0);
    N = points.ARRAYSIZE(1);
    b = 0.5;
    exps = 0.0;
    SSAValue9 = n;
    SSAValue13 = (1) <= (SSAValue9);
    SSAValue14 = (1) - (1);
    SSAValue11 = (SSAValue13) ? (SSAValue9) : (SSAValue14);
    SSAValue16 = (SSAValue11) - (1);
    SSAValue18 = (SSAValue16) + (1);
    __hpat_loop_div_7 = (SSAValue18) / (__hpat_num_pes);
    __hpat_loop_start_7 = ((__hpat_node_id) * (__hpat_loop_div_7)) + (1);
    __hpat_loop_end_7 = ((__hpat_node_id==__hpat_num_pes-1) ? SSAValue18 : (__hpat_node_id+1)*__hpat_loop_div_7);
    for ( ppip273p280 = __hpat_loop_start_7; ppip273p280 <= (int64_t)__hpat_loop_end_7; ppip273p280 += 1)
    {
        ;
        SSAValue16pp1 = X.ARRAYELEM(((ppip273p280) - (__hpat_loop_start_7)) + (1));
        parallel_ir_save_array_len_1_17 = points.ARRAYSIZE(1);
        parallel_ir_new_array_name_17_1 = j2c_array<double>::new_j2c_array_1d(NULL, parallel_ir_save_array_len_1_17);
        for ( parfor_index_1_17 = 1; parfor_index_1_17 <= (int64_t)parallel_ir_save_array_len_1_17; parfor_index_1_17 += 1)
        {
            ;
            parallel_ir_array_temp__6_18_1 = points.ARRAYELEM(parfor_index_1_17);
            SSAValue32 = (SSAValue16pp1) - (parallel_ir_array_temp__6_18_1);
            parallel_ir_array_temp__23_20_1 = SSAValue32;
            parallel_ir_array_temp_SSAValue17_22_1 = parallel_ir_array_temp__23_20_1;
            SSAValue33 = (int32_t)(2);
            pow(parallel_ir_array_temp_SSAValue17_22_1, SSAValue33);
            SSAValue34 = pow(parallel_ir_array_temp_SSAValue17_22_1, SSAValue33);
            parallel_ir_array_temp_SSAValue17_24_2 = SSAValue34;
            parallel_ir_array_temp_SSAValue18_27_1 = parallel_ir_array_temp_SSAValue17_24_2;
            SSAValue35 = -(parallel_ir_array_temp_SSAValue18_27_1);
            parallel_ir_array_temp_SSAValue18_29_2 = SSAValue35;
            parallel_ir_new_array_name_17_1.ARRAYELEM(parfor_index_1_17) = parallel_ir_array_temp_SSAValue18_29_2;
        }
        ;
        SSAValue20 = parallel_ir_new_array_name_17_1;
        SSAValue27 = (int32_t)(2);
        pow(b, SSAValue27);
        SSAValue22 = pow(b, SSAValue27);
        SSAValue28 = (double)2;
        SSAValue25pp2 = (SSAValue28) * (SSAValue22);
        parallel_ir_save_array_len_1_31 = SSAValue20.ARRAYSIZE(1);
        parallel_ir_reduction_output_35 = DBL_MAX;
        for ( parfor_index_1_31 = 1; parfor_index_1_31 <= (int64_t)parallel_ir_save_array_len_1_31; parfor_index_1_31 += 1)
        {
            ;
            parallel_ir_array_temp_SSAValue2_32_1 = SSAValue20.ARRAYELEM(parfor_index_1_31);
            SSAValue36 = (parallel_ir_array_temp_SSAValue2_32_1) / (SSAValue25pp2);
            parallel_ir_array_temp_SSAValue2_34_2 = SSAValue36;
            SSAValue20.ARRAYELEM(parfor_index_1_31) = parallel_ir_array_temp_SSAValue2_34_2;
            parallel_ir_array_temp__4_37_1 = parallel_ir_array_temp_SSAValue2_34_2;
            SSAValue37 = parallel_ir_array_temp__4_37_1;
            SSAValue38 = parallel_ir_reduction_output_35;
            SSAValue39 = (SSAValue38) < (0);
            SSAValue40 = (SSAValue37) < (0);
            SSAValue41 = !(SSAValue39);
            (SSAValue40) & (SSAValue41);
            SSAValue42 = (parallel_ir_array_temp__4_37_1) < (parallel_ir_reduction_output_35);
            SSAValue43 = (SSAValue40) & (SSAValue41);
            (SSAValue42) | (SSAValue43);
            SSAValue44 = (parallel_ir_array_temp__4_37_1) != (parallel_ir_array_temp__4_37_1);
            SSAValue45 = (parallel_ir_reduction_output_35) != (parallel_ir_reduction_output_35);
            SSAValue46 = (SSAValue42) | (SSAValue43);
            SSAValue47 = (SSAValue44) ? (parallel_ir_reduction_output_35) : (parallel_ir_array_temp__4_37_1);
            SSAValue48 = (SSAValue45) ? (parallel_ir_array_temp__4_37_1) : (parallel_ir_reduction_output_35);
            SSAValue49 = (SSAValue46) ? (SSAValue47) : (SSAValue48);
            parallel_ir_reduction_output_35 = SSAValue49;
        }
        ;
        SSAValue21 = parallel_ir_reduction_output_35;
        d = SSAValue20;
        SSAValue29 = (double)N;
        SSAValue30 = (b) * (SSAValue29);
        SSAValue23 = log(SSAValue30);;
        SSAValue24 = (SSAValue21) - (SSAValue23);
        parallel_ir_save_array_len_1_40 = d.ARRAYSIZE(1);
        parallel_ir_reduction_output_40 = DBL_MAX;
        for ( parfor_index_1_40 = 1; parfor_index_1_40 <= (int64_t)parallel_ir_save_array_len_1_40; parfor_index_1_40 += 1)
        {
            ;
            parallel_ir_array_temp__4_42_1 = d.ARRAYELEM(parfor_index_1_40);
            SSAValue50 = parallel_ir_array_temp__4_42_1;
            SSAValue51 = parallel_ir_reduction_output_40;
            SSAValue52 = (SSAValue51) < (0);
            SSAValue0 = (SSAValue50) < (0);
            SSAValue1 = !(SSAValue52);
            (SSAValue0) & (SSAValue1);
            SSAValue2 = (parallel_ir_array_temp__4_42_1) < (parallel_ir_reduction_output_40);
            SSAValue3 = (SSAValue0) & (SSAValue1);
            (SSAValue2) | (SSAValue3);
            SSAValue4 = (parallel_ir_array_temp__4_42_1) != (parallel_ir_array_temp__4_42_1);
            SSAValue5 = (parallel_ir_reduction_output_40) != (parallel_ir_reduction_output_40);
            SSAValue6 = (SSAValue2) | (SSAValue3);
            SSAValue7 = (SSAValue4) ? (parallel_ir_reduction_output_40) : (parallel_ir_array_temp__4_42_1);
            SSAValue8 = (SSAValue5) ? (parallel_ir_array_temp__4_42_1) : (parallel_ir_reduction_output_40);
            SSAValue10 = (SSAValue6) ? (SSAValue7) : (SSAValue8);
            parallel_ir_reduction_output_40 = SSAValue10;
        }
        ;
        SSAValue31pp3 = parallel_ir_reduction_output_40;
        parallel_ir_save_array_len_1_44 = d.ARRAYSIZE(1);
        parallel_ir_reduction_output_53 = 0.0;
        for ( parfor_index_1_44 = 1; parfor_index_1_44 <= (int64_t)parallel_ir_save_array_len_1_44; parfor_index_1_44 += 1)
        {
            ;
            parallel_ir_array_temp__4_45_1 = d.ARRAYELEM(parfor_index_1_44);
            SSAValue12 = (parallel_ir_array_temp__4_45_1) - (SSAValue31pp3);
            parallel_ir_array_temp__4_47_2 = SSAValue12;
            parallel_ir_array_temp_SSAValue32_49_1 = parallel_ir_array_temp__4_47_2;
            SSAValue15 = exp(parallel_ir_array_temp_SSAValue32_49_1);
            parallel_ir_array_temp_SSAValue32_51_2 = SSAValue15;
            parallel_ir_array_temp_SSAValue33_55_1 = parallel_ir_array_temp_SSAValue32_51_2;
            SSAValue17 = (parallel_ir_reduction_output_53) + (parallel_ir_array_temp_SSAValue33_55_1);
            parallel_ir_reduction_output_53 = SSAValue17;
        }
        ;
        SSAValue31 = parallel_ir_reduction_output_53;
        SSAValue25 = log(SSAValue31);;
        SSAValue26 = (SSAValue24) + (SSAValue25);
        exps = (exps) + (SSAValue26);
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
        mainFileSh << "mpiicpc -O3 -std=c++11  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << "  -lm " << std::endl;
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
        mainFileSh << "mpiicpc -O3 -std=c++11  -g   -fpic  -o " << newMainExe.str() << " " << newMain.str() << "  -lm " << std::endl;
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