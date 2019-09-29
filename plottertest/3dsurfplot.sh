#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outputdir=$scriptpath

python3 $scriptpath/../plotter.py \
    --csvlog "$scriptpath/3dsurfplottest.csv" \
    --graphpath "$outputdir/3dsurfplottest.pdf" \
    --show

