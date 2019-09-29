#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outputdir=$scriptpath

python3 $scriptpath/../plotter.py \
    --csvlog "$scriptpath/3dsurfplottest.csv" \
    --graphpath "$outputdir/3dsurfplottest.pdf" \
    --only-tags "label-a" \
    --title "3dsurf Test" \
    --xlabel "X-Label" \
    --ylabel "Y-Label" \
    --custom-legend-labels "custom-label-a" \
    --color-theme "dracula" \
    --color-offset 1 \
    --show

