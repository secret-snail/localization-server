#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
logdir=$scriptpath/../results/
plotdir=$scriptpath/../plots/

#SHOW="--show"
SHOW=""

# relies on:
#   emp_float_benchmark_run.sh
#   emp_float_benchmark_dataobl_run.sh

sed -i 's/grouped_bar,operator(),emp_float_gn_time_vs_point/xy,operator(),emp_float_gn_time_vs_point/g' $scriptpath/../results/gn_emp_float_eth3d_bench_short.log
sed -i 's/grouped_bar,operator(),emp_float_lm_time_vs_point/xy,operator(),emp_float_lm_time_vs_point/g' $scriptpath/../results/lm_emp_float_eth3d_bench_short.log
sed -i 's/grouped_bar,operator(),emp_float_gn_time_vs_point_dataobl/xy,operator(),emp_float_gn_time_vs_point_dataobl/g' $scriptpath/../results/gn_emp_float_eth3d_bench_dataobl_short.log
sed -i 's/grouped_bar,operator(),emp_float_lm_time_vs_point_dataobl/xy,operator(),emp_float_lm_time_vs_point_dataobl/g' $scriptpath/../results/lm_emp_float_eth3d_bench_dataobl_short.log

python3 $scriptpath/../scripts/plotter-do-vs-sil.py \
     --csvlog \
        "$logdir/gn_emp_float_eth3d_bench_short.log" \
        "$logdir/lm_emp_float_eth3d_bench_short.log" \
        "$logdir/gn_emp_float_eth3d_bench_dataobl_short.log" \
        "$logdir/lm_emp_float_eth3d_bench_dataobl_short.log" \
     --graphpath "$plotdir/loopleak_vs_dataobl.pdf" \
     --title "Comparing Data Oblivious and
Single Iteration Localization Runtime" \
     --only-tags \
        "emp_float_gn_time_vs_points_dataobl" \
        "emp_float_lm_time_vs_points_dataobl" \
        "emp_float_gn_time_vs_points" \
        "emp_float_lm_time_vs_points" \
     --custom-legend-labels \
        "GN DO" \
        "LM DO" \
        "GN SIL" \
        "LM SIL" \
     --xlabel "Number of Features" \
     --ylabel "Runtime (s)" \
     --fig-w 4 \
     --fig-h 3 \
     --log-scale \
     --color-theme "dracula" \
     $SHOW
