#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

NUM_FRAMES=1
NUM_TRIALS=3

# long run only (256 points)
$scriptpath/../build/bin/emp_float_eth3d_bench gn $NUM_FRAMES $NUM_TRIALS 256 | tee $scriptpath/../results/gn_emp_float_eth3d_bench_long_latency.log
$scriptpath/../build/bin/emp_float_eth3d_bench lm $NUM_FRAMES $NUM_TRIALS 256 | tee $scriptpath/../results/lm_emp_float_eth3d_bench_long_latency.log

sed -i 's/emp_float_gn_time_vs_points_per_loc_itr/emp_float_gn_time_vs_points_per_loc_itr_latency/g' $scriptpath/../results/gn_emp_float_eth3d_bench_long_latency.log
sed -i 's/emp_float_lm_time_vs_points_per_loc_itr/emp_float_lm_time_vs_points_per_loc_itr_latency/g' $scriptpath/../results/lm_emp_float_eth3d_bench_long_latency.log
