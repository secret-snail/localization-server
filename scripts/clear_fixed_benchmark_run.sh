#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

killall gn_fixed_server_benchmark
killall lm_fixed_server_benchmark

$scriptpath/../build/fixed_server_benchmark 0 1 6 > $scriptpath/../results/fixed_server_benchmark_short.log &
$scriptpath/../build/fixed_server_benchmark 1 1 6
