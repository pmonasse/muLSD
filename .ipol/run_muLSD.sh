#!/bin/bash
# Arguments: scales (>=0)
# -zoom: zoom factor
# -color: handle color image?

scales=$1

$bin/build/muLSD -s $scales $input_0 mulsd.txt &&
$bin/build2/draw_lines mulsd.txt lines.png &&
(! [ -x composite ] || composite lines.png $input_0 linesOver.bmp)
