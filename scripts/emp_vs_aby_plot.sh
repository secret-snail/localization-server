#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
logdir=$scriptpath/../results/
plotdir=$scriptpath/../plots/

#SHOW="--show"
SHOW=""

# relies on:
#   aby_float_benchmark_run.sh
#   emp_float_benchmark_run.sh

sed -i 's/xy,operator(),emp_float_gn_time_vs_point/grouped_bar,operator(),emp_float_gn_time_vs_point/g' $scriptpath/../results/gn_emp_float_eth3d_bench_short.log
sed -i 's/xy,operator(),emp_float_lm_time_vs_point/grouped_bar,operator(),emp_float_lm_time_vs_point/g' $scriptpath/../results/lm_emp_float_eth3d_bench_short.log
sed -i 's/xy,operator(),aby_bool_float_gn_time_vs_point/grouped_bar,operator(),aby_bool_float_gn_time_vs_point/g' $scriptpath/../results/gn_aby_float_eth3d_bench_short.log
sed -i 's/xy,operator(),aby_yao_float_gn_time_vs_point/grouped_bar,operator(),aby_yao_float_gn_time_vs_point/g' $scriptpath/../results/gn_aby_float_eth3d_bench_short.log
sed -i 's/xy,operator(),aby_bool_float_lm_time_vs_point/grouped_bar,operator(),aby_bool_float_lm_time_vs_point/g' $scriptpath/../results/lm_aby_float_eth3d_bench_short.log
sed -i 's/xy,operator(),aby_yao_float_lm_time_vs_point/grouped_bar,operator(),aby_yao_float_lm_time_vs_point/g' $scriptpath/../results/lm_aby_float_eth3d_bench_short.log

# may need to change plot tags from xy to grouped_bar
# comment out twelve points :%s/\(.*\),12,/#\1,12,/g
python3 $scriptpath/../scripts/plotter-emp-vs-aby.py \
     --csvlog \
        "$logdir/gn_aby_float_eth3d_bench_short.log" \
        "$logdir/lm_aby_float_eth3d_bench_short.log" \
        "$logdir/gn_emp_float_eth3d_bench_short.log" \
        "$logdir/lm_emp_float_eth3d_bench_short.log" \
     --graphpath "$plotdir/emp_vs_aby.pdf" \
     --title "ABY and EMP Localization Runtime" \
     --only-tags \
        "aby_yao_float_gn_time_vs_points" \
        "aby_bool_float_gn_time_vs_points" \
        "aby_yao_float_lm_time_vs_points" \
        "aby_bool_float_lm_time_vs_points" \
        "emp_float_gn_time_vs_points" \
        "emp_float_lm_time_vs_points" \
     --custom-legend-labels \
        "ABY Yao GN" \
        "ABY Bool GN" \
        "ABY Yao LM" \
        "ABY Bool LM" \
        "EMP GN" \
        "EMP LM" \
     --xlabel "Number of Features" \
     --ylabel "Runtime (s)" \
     --fig-w 4 \
     --fig-h 5 \
     --log-scale \
     --horizontal \
     --color-theme "dracula" \
     $SHOW
