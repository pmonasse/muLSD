#!/bin/bash
# Arguments: scales (>=0)
# -zoom: zoom factor
# -color: handle color image?

scales=$1

$bin/build/muLSD -s $scales $input_0 mulsd.txt
