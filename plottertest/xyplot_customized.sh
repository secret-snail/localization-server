#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outputdir=$scriptpath

python3 $scriptpath/../plotter.py \
    --csvlog "$scriptpath/xyplottest.csv" \
    --graphpath "$outputdir/xyplottest.pdf" \
    --only-tags "label-a" "label-b" \
    --title "xy Test" \
    --xlabel "X-Label" \
    --ylabel "Y-Label" \
    --custom-legend-labels "label-a" "label-b" \
    --color-theme "dracula" \
    --color-offset 1 \
    --show

