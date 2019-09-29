#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outputdir=$scriptpath

python3 $scriptpath/../plotter.py \
    --csvlog "$scriptpath/3dbarplottest.csv" \
    --graphpath "$outputdir/3dbarplottest.pdf" \
    --show

