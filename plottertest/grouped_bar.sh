#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outputdir=$scriptpath

python3 $scriptpath/../plotter.py \
    --csvlog "$scriptpath/grouped_bar.csv" \
    --graphpath "$outputdir/grouped_bar.pdf" \
    --show

