# Secure Localization Server
[![devcontainer](https://github.com/secret-snail/localization-server/actions/workflows/devcontainer.yml/badge.svg)](https://github.com/secret-snail/localization-server/actions/workflows/devcontainer.yml)
[![container](https://github.com/secret-snail/localization-server/actions/workflows/docker-image.yml/badge.svg)](https://github.com/secret-snail/localization-server/actions/workflows/docker-image.yml)

Privacy preserving localization based on secure multiparty computation (MPC).

## Build via Container
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
ABY Dependencies:
`sudo apt install g++ make cmake libgmp-dev libssl-dev libboost-all-dev`

EMP Dependencies:
`sudo apt install software-properties-common cmake git build-essential libssl-dev`

Build and test the code:
```
git submodule update --init --recursive
mkdir build && cd build
cmake ..
make
ctest
```

ABY testing requires -DCMAKE_BUILD_TYPE=Release, tests fail when the sanitizers
are turned on because there appears to be memory leaks in the ABY library.

Note - circuits cannot be built on the fly. must be fully specified then executed.
This means if control flow requires some secret data, circuit must be broken and
intermediate ciphertext stored as secret share.
[Developer Guide](https://www.informatik.tu-darmstadt.de/media/encrypto/encrypto_code/abydevguide.pdf)
[Reusing Computation](https://github.com/encryptogroup/ABY/issues/167)

## Experimental Evaluation
First, download the
[eth3d dataset](https://www.eth3d.net/datasets#high-res-multi-view) with the
following commands. These will specifically download the high-res multi-view
undistorted images and ground truth scan evaluation from eth3d.
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
docker/build.sh

mkdir -p results

docker run -it --rm --init \
  --net=host \
  --name snail-tester \
  --volume "$(pwd)/results":/snail/results \
  snail-server bash
```

Run the experiments. Note some of the commands must be run outside the container
due to network setup requirements. They are marked with `# outside container`.

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

Plot the results.

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
Make sure the snail can see the marker in the first 10 seconds, otherwise the
snail will need to be restarted.


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
