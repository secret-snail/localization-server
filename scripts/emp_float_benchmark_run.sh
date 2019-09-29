#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

NUM_FRAMES=1
NUM_TRIALS=3

# short run (12 points)
$scriptpath/../build/bin/emp_float_eth3d_bench gn $NUM_FRAMES $NUM_TRIALS 12 | tee $scriptpath/../results/gn_emp_float_eth3d_bench_short.log
$scriptpath/../build/bin/emp_float_eth3d_bench lm $NUM_FRAMES $NUM_TRIALS 12 | tee $scriptpath/../results/lm_emp_float_eth3d_bench_short.log

# large run (256 points)
$scriptpath/../build/bin/emp_float_eth3d_bench gn $NUM_FRAMES $NUM_TRIALS 256 | tee $scriptpath/../results/gn_emp_float_eth3d_bench_long.log
$scriptpath/../build/bin/emp_float_eth3d_bench lm $NUM_FRAMES $NUM_TRIALS 256 | tee $scriptpath/../results/lm_emp_float_eth3d_bench_long.log