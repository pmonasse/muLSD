# Procedure for benchmarking with York Urban Line Segment Database

## Preparation (only once)

Steps 2 and 3 are optional, the ground truth txt files are already present

1. Download the dataset:
https://www.elderlab.yorku.ca/?smd_process_download=1&download_id=8288
2. Ground truth is stored in mat files, we want txt. To convert, go to folder containing the data  and run:
```shell
matlab -batch "export_GT; exit"
```
This creates the files P???????_gt.txt
3. Copy these files to the current folder (benchmark)
```shell
mv /path/P*_gt.txt .
```
4. Link the image files here:
```shell
ln -s $(find /path -name \*.pgm) .
```

## Benchmarking

1. Standard run with muLSD:
```shell
make clean
make -s -j 2
```
(replace 2 above with higer number if you have more cores to put to work).
It displays the precision, recall and IOU.
2. Run with non-default parameters:
```shell
make clean
LSD="../build/muLSD -s 1 -g 2" make -s -j 2
```
3. Run with baseline LSD: change the EXT variable because LSD reads PGM files:
```shell
make clean
LSD="../lsd_1.6/lsd" EXT=pgm make -s -j 2
```