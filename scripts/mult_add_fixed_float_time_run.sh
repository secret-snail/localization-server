#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

killall emp_float_vs_fixed_benchmark_*
$scriptpath/../build/bin/emp_float_vs_fixed_benchmark_add 0 > $scriptpath/../results/emp_float_vs_fixed_benchmark_add.log &
$scriptpath/../build/bin/emp_float_vs_fixed_benchmark_add 1

killall emp_float_vs_fixed_benchmark_*
$scriptpath/../build/bin/emp_float_vs_fixed_benchmark_mul 0 > $scriptpath/../results/emp_float_vs_fixed_benchmark_mul.log &
$scriptpath/../build/bin/emp_float_vs_fixed_benchmark_mul 1
