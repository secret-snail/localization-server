# Scripts for Benchmarking

This directory contains scripts for benchmarking the performance of the
localization libraries and can be used to reproduce the results in the paper. At
a high level, each experiment has a `*_run.sh` script which runs one or more
benchmarks, piping the console output to files in `results/`. The respective
`*_plot.sh` script reads the results and generates plots by calling the
appropriate python script in this directory and stores the output in `plots/`.

The following mapping shows the relationship between the scripts and the
logical experiments performed:

- Experiment 1: The proposed approach -- single iteration localization.
  - `scripts/emp_float_benchmark_run.sh`
  - `scripts/aby_float_benchmark_run.sh`
  - `scripts/mult_add_fixed_float_time_run.sh`
- Experiment 2: A baseline data oblivious adaptation of localization.
  - `scripts/emp_float_benchmark_dataobl_run.sh`
- Experiment 3: Single iteration localization where 5ms of network latency has
been introduced between the two servers.
  - `scripts/emp_float_benchmark_run_latency.sh`

The following mapping shows the relationship between the experiments and figures
in the paper:

- Experiment 1 generates Figures 2, 3, 4, 5, and 6.
- Experiment 2 generates Figures 6.
- Experiment 3 generates Figures 7a and 7b.

See the top level [README.md](../README.md) for instructions on running the
experiments.
