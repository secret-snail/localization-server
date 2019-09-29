#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
logdir=$scriptpath/../results/
plotdir=$scriptpath/../plots/

#SHOW="--show"
SHOW=""


# Plot error
python3 $scriptpath/../scripts/plotter.py \
     --csvlog "$logdir/gn_fixed_eth3d_bench_short.log" \
     --graphpath "$plotdir/gn_fixed_eth3d_bench_error.pdf" \
     --title "GN Error vs Ground Truth (plaintext)" \
     --only-tags "opencv_error_vs_numpts" \
        "gn_float_error_vs_numpts" \
        "gn_baset_error_vs_numpts" \
     --custom-legend-labels "opencv (float)" \
        "float" \
        "fixed (64bit)" \
     --xlabel "Number of Points" \
     --ylabel "Error (2norm)" \
     --color-theme "dracula" \
     --fig-w 5 \
     --fig-h 4 \
     $SHOW

# Plot number of iterations until convergence
python3 $scriptpath/../scripts/plotter.py \
     --csvlog "$logdir/gn_fixed_eth3d_bench_short.log" \
     --graphpath "$plotdir/gn_fixed_eth3d_bench_iterations.pdf" \
     --title "GN Iterations (plaintext)" \
     --only-tags \
        "gn_float_iterations_vs_numpts" \
        "gn_baset_iterations_vs_numpts" \
     --custom-legend-labels \
        "float" \
        "fixed (64bit)" \
     --xlabel "Number of Points" \
     --ylabel "Iterations" \
     --color-theme "dracula" \
     --fig-w 5 \
     --fig-h 4 \
     $SHOW
