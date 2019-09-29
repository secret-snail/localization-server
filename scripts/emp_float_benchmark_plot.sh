#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
logdir=$scriptpath/../results/
plotdir=$scriptpath/../plots/

#SHOW="--show"
SHOW=""

# relies on:
#   emp_float_benchmark_run.sh
#   emp_float_benchmark_run_latency.sh

sed -i 's/,emp_float_gn_time_vs_points_per_loc_itr,/,emp_float_gn_time_vs_points_per_loc_itr_latency,/g' $scriptpath/../results/gn_emp_float_eth3d_bench_long_latency.log
sed -i 's/,emp_float_lm_time_vs_points_per_loc_itr,/,emp_float_lm_time_vs_points_per_loc_itr_latency,/g' $scriptpath/../results/lm_emp_float_eth3d_bench_long_latency.log

python3 $scriptpath/../scripts/insetplotter.py \
     --csvlog \
        "$logdir/gn_emp_float_eth3d_bench_long.log" \
        "$logdir/lm_emp_float_eth3d_bench_long.log" \
        "$logdir/gn_emp_float_eth3d_bench_long_latency.log" \
        "$logdir/lm_emp_float_eth3d_bench_long_latency.log" \
     --graphpath "$plotdir/emp_float_runtime_long.pdf" \
     --title "Comparing Optimization
Algorithms at Large Input Sizes" \
     --only-tags \
        "emp_float_gn_time_vs_points_per_loc_itr" \
        "emp_float_lm_time_vs_points_per_loc_itr" \
        "emp_float_gn_time_vs_points_per_loc_itr_latency" \
        "emp_float_lm_time_vs_points_per_loc_itr_latency" \
     --custom-legend-labels \
        "GN (0 ms)" \
        "LM (0 ms)" \
        "GN (5 ms)" \
        "LM (5 ms)" \
     --xlabel "Number of Features" \
     --ylabel "Runtime per Iteration (s)" \
     --fig-w 4 \
     --fig-h 3 \
     --color-theme "dracula" \
     $SHOW
