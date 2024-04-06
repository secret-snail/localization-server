# Secure Localization Server
[![devcontainer](https://github.com/secret-snail/localization-server/actions/workflows/devcontainer.yml/badge.svg)](https://github.com/secret-snail/localization-server/actions/workflows/devcontainer.yml)
[![container](https://github.com/secret-snail/localization-server/actions/workflows/docker-image.yml/badge.svg)](https://github.com/secret-snail/localization-server/actions/workflows/docker-image.yml)

**Privacy preserving localization using secure multiparty computation.**

Visual localization is a computer vision task by which the position and
orientation of a camera is determined from an image and environmental map. This
project implements a secure localization server that uses MPC to compute a
camera pose without revealing the image or map. For more information, see the
paper on [arxiv](https://arxiv.org/abs/2403.14916), published at
[PETs 2023](https://petsymposium.org/2024/paperlist.php).

## Repository Structure

This repository stores external dependencies as submodules in the `extern`
directory. The core localization libraries for each MPC protocol are stored in
the `src` directory. The `test` directory contains unit tests and benchmarks for
the localization libraries. The `.vscode` and `.devcontainer` directories
automatically configure a working development environment and continuous
deployment is set up with GitHub Actions in `.github`.

## Build via Container

Building the container requires compiling OpenCV, which may take 30+
minutes.

```bash
docker/build.sh
```

Run an interactive shell inside the conatiner:

```bash
docker run -it --rm --init \
  --net=host \
  --name snail-server \
  snail-server bash
```

## Build From Source

If not using the provided devcontainer, install the following dependencies:

```bash
sudo apt install g++ make cmake libgmp-dev libssl-dev libboost-all-dev  # ABY
sudo apt install software-properties-common cmake git build-essential libssl-dev  # EMP
```

Build and test the code:

```bash
git submodule update --init --recursive
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
ctest
```

Special notes about ABY, one of the MPC libraries used in this project:
- ABY tests and experiments are memory intensive, requiring up to 19 GB of
memory (resident set size).
- ABY tests require `-DCMAKE_BUILD_TYPE=Release`. The tests fail when the
sanitizers are turned on because there appear to be memory leaks in the external
ABY library. Due to this, tracking memory leaks in the ABY localization
library in *this* repo is difficult thus likely also leaks, though an effort was
made to avoid this. This is not true for the EMP localization library.
- ABY circuits cannot be built on the fly. They must be fully specified then
executed. This means if control flow requires some secret data, circuit must be
broken and intermediate ciphertext stored as secret share.
[Developer Guide](https://www.informatik.tu-darmstadt.de/media/encrypto/encrypto_code/abydevguide.pdf)
[Reusing Computation](https://github.com/encryptogroup/ABY/issues/167)

Special notes about EMP, the other MPC library used in this project:
- A machine with AES hardware acceleration (for example AES-NI on Intel
hardware) is required for performance reasons.

### External Dependencies

External source dependencies in `extern`:

- [ABY](https://github.com/encryptogroup/ABY)
- [Catch2](https://github.com/catchorg/Catch2)
- [EMP-Toolkit](https://github.com/emp-toolkit)
- [OpenCV](https://github.com/opencv/opencv)
- [Fixed Point Math Library](https://github.com/MikeLankamp/fpm) (direct import at `src/common/fixed_point.h`)
- [Singular Value Decomposition](https://numerical.recipes/) (direct import at `test/common-test/cleartext-ref/svd.hpp`)

Other dependencies:

- docker
- cmake
- gcc
- clang-format
- python3
- matplotlib

## Reproduce Results

First, download the
[eth3d dataset](https://www.eth3d.net/datasets#high-res-multi-view) with the
following commands on the host (outside the container). These will specifically
download the high-res multi-view undistorted images and ground truth scan
evaluation from eth3d. Note this requires ~20GB of storage.

```bash
sudo apt-get install p7zip-full
mkdir ./data-eth3d
curl https://www.eth3d.net/data/multi_view_training_dslr_undistorted.7z -o im.7z
7z x im.7z -o./data-eth3d/
rm im.7z
curl https://www.eth3d.net/data/multi_view_training_dslr_scan_eval.7z -o gt.7z
7z x gt.7z -o./data-eth3d/
rm gt.7z
```

The easiest way to run the experiments is via the provided Docker container.
Build and start the container with the following commands.

```bash
docker/build.sh  # 30 minutes

mkdir -p results

docker run -it --rm --init \
  --net=host \
  --name snail-tester \
  --volume "$(pwd)/data-eth3d":/snail/data-eth3d \
  --volume "$(pwd)/results":/snail/results \
  snail-server bash
```

Run `ctest` inside the `build` directory of the container to verify the build.
The tests should pass in roughly 15 minutes. Machine with less than 19 GB of
memory may fail some of the ABY tests due to memory constraints. Machines with
sufficient memory may still occasionally fail the `ABY large SVD` test due to
memory errata.

Run the experiments with the following commands. Note some of the commands must
be run outside the container due to network setup requirements. They are marked
with `# outside container`. Each experiment outputs information to the console
for example the image being tested, how many points are being used, the ground
truth pose and the pose as estimated by OpenCV vs. the MPC protocol. The
console output is piped to a file in `results` and lines prefixed by `SeNtInAl`
are used to parse the results in subsequent plotting. The characters which look
like `...:::555/44/33/2/1` correspond to the different phases and iterations of
the SVD algorithm. Lines which look like "Trial x with y randomly selected
points..." show the progress of the experiment.


```bash
# ABY and EMP single-iteration localization tests.
sed -i 's/#define PPL_FLOW .*/#define PPL_FLOW PPL_FLOW_SiSL/' src/common/privacyconf.h
(cd build/ && make)
LAT=0msec source scripts/network_setup.sh # outside container
scripts/emp_float_benchmark_run.sh # 2.5 hours
scripts/aby_float_benchmark_run.sh # 7.5 hours
scripts/mult_add_fixed_float_time_run.sh # 30 seconds

# EMP data-oblivious tests.
sed -i 's/#define PPL_FLOW .*/#define PPL_FLOW PPL_FLOW_DO/' src/common/privacyconf.h
(cd build/ && make)
LAT=0msec source scripts/network_setup.sh # outside container
scripts/emp_float_benchmark_dataobl_run.sh # 38 hours

# EMP single-iteration localization tests with network latency.
sed -i 's/#define PPL_FLOW .*/#define PPL_FLOW PPL_FLOW_SiSL/' src/common/privacyconf.h
(cd build/ && make)
LAT=5msec source scripts/network_setup.sh # outside container
scripts/emp_float_benchmark_run_latency.sh # 3 hours

source scripts/network_teardown.sh # outside container
exit # stop the container
```

Install dependencies for the plotter scripts (requires python3 and pip3).

```bash
./plotter_deps.sh
```

Plot the results. The outputs are stored in `plots/` and can be matched to
figures in the paper (see [claims](claims.md)). The time required to exactly
reproduce the experiments from the paper is large thus the number of averaged
trials for each experiment has been reduced from results shown in the paper.
To run an exact reproduction, the scripts in `scripts/*_run.sh` may be edited to
run 5 trials and 8 frames.

```bash
mkdir plots
scripts/emp_vs_aby_plot.sh
scripts/loopleak_vs_dataobl_plot.sh
scripts/emp_float_benchmark_plot.sh
scripts/netio_plot.sh
scripts/num_arith_ops_plot.sh
scripts/mult_add_fixed_float_time_plot.sh
```

## Robotic Snail Demo

*Requires raspberry pi, see [snail repo](https://github.com/secret-snail/snail)*

Currently, the demo requires `privacy_conf.h` to have this set:
```
#define PPL_FLOW PPL_FLOW_LOOP_LEAK
```

To run the snail demo, first `docker/build.sh` the container. Then run
`docker/alice.sh` and `docker/bob.sh` in two terminals on the server(s).
On the [snail](https://github.com/secret-snail/snail), run
`sudo ./build/bin/visp_snail --secure`.

## Running on Separate Machines

Set the machines IP addresses in `test/emp-float/client.cpp` and
`src/emp-float-server/server-lm.cpp` replacing the localhost IP.

```bash
cd build/bin
./emp_float_client
```

```bash
cd build/bin
./lm_emp_float_server 0
```

```bash
cd build/bin
./lm_emp_float_server 1
```
