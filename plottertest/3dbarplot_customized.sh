#!/bin/bash
scriptpath="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outputdir=$scriptpath

python3 $scriptpath/../plotter.py \
    --csvlog "$scriptpath/3dbarplottest.csv" \
    --graphpath "$outputdir/3dbarplottest.pdf" \
    --only-tags "label-a" "label-b" "label-c" \
    --title "3dbar Test" \
    --xlabel "X-Label" \
    --ylabel "Y-Label" \
    --custom-legend-labels "custom-label-a" "custom-label-b" "custom-label-c" \
    --color-theme "dracula" \
    --color-offset 1 \
    --show

