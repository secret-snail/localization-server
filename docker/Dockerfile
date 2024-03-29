FROM ubuntu:22.04

RUN apt update && apt upgrade -y && DEBIAN_FRONTEND=noninteractive apt install -y sudo
RUN DEBIAN_FRONTEND=noninteractive apt install -y g++ make cmake libgmp-dev libssl-dev libboost-all-dev

RUN DEBIAN_FRONTEND=noninteractive apt install -y software-properties-common cmake git build-essential libssl-dev xxd

WORKDIR /snail

# build opencv first for caching
COPY ./extern/opencv ./extern/opencv
COPY ./extern/opencv_contrib ./extern/opencv_contrib
RUN cd extern/opencv && \
  mkdir build && \
  cmake -S ./ -B ./build/ && \
  cmake --build ./build/ --parallel 8 && \
  cd build && sudo make install

COPY ./ ./

RUN mkdir build && \
  cd build && \
  cmake -DCMAKE_BUILD_TYPE=Release .. && \
  make -j 8
