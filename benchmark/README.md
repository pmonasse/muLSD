# Procedure for benchmarking with York Urban Line Segment Database

Notice we rely on YorkUrban (1)[^1] with enriched annotations (2)[^2].

## Preparation (only once)

Steps 2 and 3 are optional, the ground truth txt files are already present

1. Download the dataset:
https://www.dropbox.com/scl/fi/lzm76drgp97mmy0fpygz7/YorkUrban-LineSegment.zip?e=1&rlkey=emo2b65onipofxikrdzd2kufm
2. Ground truth is stored in mat files, we want txt. To convert, go to folder LineSegmentAnnotation and run:
```shell
cd YorkUrban-LineSegment/LineSegmentAnnotation
matlab -batch "oldpath=path; path(oldpath,'../..'); export_GT; exit"
cd ../..
```
This creates the files P???????_gt.txt

3. Copy these files to the current folder (benchmark)
```shell
mv YorkUrban-LineSegment/LineSegmentAnnotation/P*_gt.txt .
```
4. Link the image files here:
```shell
ln -s $(find YorkUrban-LineSegment/P* -name \*.jpg) .
```

## Benchmarking

1. Standard run with muLSD:
```shell
make clean
make -s -j 2
```
(replace 2 above with higer number if you have more CPU cores to put to work).
It displays the precision, recall and IOU.

2. Run with non-default parameters:
```shell
make clean
LSD="../Build/muLSD -s 1 -g 2" make -s -j 2
```
3. Run with baseline LSD: change the EXT variable because LSD reads PGM files:
```shell
make clean
LSD="../lsd_1.6/lsd" EXT=pgm make -s -j 2
```

## References

[^1] Original YorkUrban dataset
```
@Inbook{Denis2008,
  author="Denis, P. and Elder, J.H. and Estrada, F.J.",
  title="Efficient Edge-Based Methods for Estimating {M}anhattan Frames in Urban Imagery",
  bookTitle="10th European Conference on Computer Vision",
  year="2008",
  pages="197--210",
  doi="10.1007/978-3-540-88688-4_15"
}
```

[^2] Augmented YorkUrban dataset
```
@article{Namgyu2017TPAMI, 
  author={Cho, N.-G. and Yuille, A. and Lee, S.-W.}, 
  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence}, 
  title={A Novel Linelet-based Representation for Line Segment Detection}, 
  year={2018}, 
  volume={40}, 
  number={5}, 
  pages={1195--1208}, 
  doi={10.1109/TPAMI.2017.2703841}
}
```
