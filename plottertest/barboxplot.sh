#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outputdir=$scriptpath

python3 $scriptpath/../plotter.py \
    --csvlog "$scriptpath/barboxplottest.csv" \
    --graphpath "$outputdir/barboxplottest.pdf" \
    --show

