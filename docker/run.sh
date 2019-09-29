#!/bin/bash

scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $scriptpath

cmd="./build/bin/lm_emp_float_server 0 &
./build/bin/lm_emp_float_server 1"

docker run -it --rm \
  --init \
  -p 8097:8097 \
  -p 8114:8114 \
  --name snail-server \
  snail-server bash -c "$cmd"
