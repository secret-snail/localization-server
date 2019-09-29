#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
logdir=$scriptpath/../results/
plotdir=$scriptpath/../plots/

#SHOW="--show"
SHOW=""

# relies on:
#   emp_float_benchmark_run.sh

# NOTE: For GB you need to modify plotter script

python3 $scriptpath/../scripts/plotter-netio.py \
     --csvlog \
        "$logdir/gn_emp_float_eth3d_bench_long.log" \
        "$logdir/lm_emp_float_eth3d_bench_long.log" \
     --graphpath "$plotdir/netio.pdf" \
     --title "Comparing Network IO
at Large Input Sizes" \
     --only-tags \
        "gn_bytes_rx_per_itr" \
        "lm_bytes_rx_per_itr" \
        "gn_bytes_rx_per_itr_per_feat" \
        "lm_bytes_rx_per_itr_per_feat" \
     --custom-legend-labels \
        "GN" \
        "LM" \
        "GN per Feature" \
        "LM per Feature" \
     --xlabel "Number of Features" \
     --ylabel "Bytes per Iteration" \
     --fig-w 4 \
     --fig-h 3 \
     --color-theme 'dracula' \
     $SHOW
