#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
logdir=$scriptpath/../results/
plotdir=$scriptpath/../plots/

#SHOW="--show"
SHOW=""


# To plot the add multiply together
python3 $scriptpath/../scripts/plotter-mult-add.py \
    --csvlog \
       "$logdir/emp_float_vs_fixed_benchmark_add.log" \
       "$logdir/emp_float_vs_fixed_benchmark_mul.log" \
    --graphpath "$plotdir/emp_float_vs_fixed_benchmark_add_mul.pdf" \
    --title "Comparing Arithmetic Operations
Across Data Representations" \
    --only-tags \
        "EMP_ADD" \
        "EMP_MUL" \
    --custom-legend-labels \
        "Add" \
        "Multiply" \
    --xlabel "" \
    --ylabel "" \
    --fig-w 4 \
    --fig-h 3 \
    --rotate_x_labels \
    --color-theme "dracula" \
    $SHOW
    #--log-scale \
    #--only-tags "opencv_error_vs_numpts" \
    #--custom-legend-labels "opencv (float)" \
    #   "float" \
    #   "fixed (64bit)" \



# To plot the add multiply graphs separately
#python3 $scriptpath/../scripts/plotter.py \
#    --csvlog \
#       "$logdir/emp_float_vs_fixed_benchmark_add.log" \
#    --graphpath "$plotdir/emp_float_vs_fixed_benchmark_add.pdf" \
#    --title "Float vs Fixed Operation Speed (10k Operations)" \
#    --xlabel "" \
#    --ylabel "Time (s)" \
#    --color-theme "dracula" \
#    --fig-w 5 \
#    --fig-h 2.5 \
#    --rotate_x_labels \
#    $SHOW
#    #--only-tags "opencv_error_vs_numpts" \
#    #--custom-legend-labels "opencv (float)" \
#    #   "float" \
#    #   "fixed (64bit)" \
#
#python3 $scriptpath/../scripts/plotter.py \
#    --csvlog \
#       "$logdir/emp_float_vs_fixed_benchmark_mul.log" \
#    --graphpath "$plotdir/emp_float_vs_fixed_benchmark_mul.pdf" \
#    --xlabel "" \
#    --ylabel "Time (s)" \
#    --color-theme "dracula" \
#    --fig-w 5 \
#    --fig-h 2.5 \
#    --rotate_x_labels \
#    --color-offset 1 \
#    $SHOW
#    #--title "Float vs Fixed Operation Speed" \
