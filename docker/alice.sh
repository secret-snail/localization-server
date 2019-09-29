#!/bin/bash

scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $scriptpath

cmd="./build/bin/lm_emp_float_server 0"

# need --init so ctrl-c signal is passed to proc correctly
# need --net=host if two containers on same machine talking
#     through localhost, because you can't open same port (8080)
#     for two containers.
#     otherwise open the following ports:
    #-p 8080:8080 \
    #-p 8097:8097 \

docker run -it --rm --init \
  --net=host \
  --name snail-server-alice \
  snail-server $cmd
