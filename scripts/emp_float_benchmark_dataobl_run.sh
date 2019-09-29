#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# only runs the data oblivious implementation.
# dont forget to set privacyconf.h PPL_FLOW to data oblivious.

# short run
$scriptpath/../build/bin/emp_float_eth3d_bench gn 3 2 12 | tee $scriptpath/../results/gn_emp_float_eth3d_bench_dataobl_short.log
$scriptpath/../build/bin/emp_float_eth3d_bench lm 3 2 12 | tee $scriptpath/../results/lm_emp_float_eth3d_bench_dataobl_short.log

