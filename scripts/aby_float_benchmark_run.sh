#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

NUM_FRAMES=1
NUM_TRIALS=3

$scriptpath/../build/bin/aby_float_eth3d_bench gn $NUM_FRAMES $NUM_TRIALS 12 | tee $scriptpath/../results/gn_aby_float_eth3d_bench_short.log
$scriptpath/../build/bin/aby_float_eth3d_bench lm $NUM_FRAMES $NUM_TRIALS 12 | tee $scriptpath/../results/lm_aby_float_eth3d_bench_short.log