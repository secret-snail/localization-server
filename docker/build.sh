#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd $scriptpath/..

docker build -f ./docker/Dockerfile --tag=snail-server .
