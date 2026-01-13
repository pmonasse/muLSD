## muLSD (a Multiscale Line Segment Detector done right)

Version 0.9, 01/13/2026.

Future versions: <https://github.com/pmonasse/muLSD>

### Authors
Pascal Monasse <pascal.monasse@enpc.fr>,
Rahima Djahel <rdjahel@gmail.com>,
Yohann Salaun <yohann.salaun@imagine.enpc.fr>

### Introduction
This code implements a multi-scale line segment detector, called muLSD. This is an extension of LSD[^1] yielding better results when applied to large images. Such an extension, called MLSD, was already performed in [^2], and implemented in [^3]. However, this was flawed and the current one is based on the same ideas but fixes several problems.

### Building
*Requirements:*

  - CMake (available at https://cmake.org/download/)
  - C++ compiler

*Build instructions:*

- Unix, MacOS:
  ```
  $ cd /path_to_this_file/
  $ cmake -DCMAKE_BUILD_TYPE:bool=Release -S . -B Build
  $ cmake --build Build
  ```
- Windows with MinGW:
  ```
  $ cd /path_to_this_file/
  $ cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE:bool=Release S . -B Build
  $ cmake --build Build
  ```

### Usage
```
Build/muLSD [options] imgIn.png out.txt
-s, --scales=ARG nb scales (0=automatic) (0)
-g, --gradient=ARG Min gradient norm (0=automatic) (0)
```
Options:

- -s: 1=single-scale, 2=single scale and half-scale, etc. 0=automatic.
- -g: 0 means the threshold is based on standard deviation of gradient norm.

Example:
```
Build/muLSD data/grid.png grid.txt
```
Generated file `grid.txt` must coincide with reference `data/grid_lines.txt`.

### Files
cluster.cpp(\*)     image.cpp  lsd.c          mulsd.hpp(\*)
cluster.hpp(\*)     image.h    lsd.h          segment.cpp(\*)
cmdLine.h          io_png.c    main.cpp(\*)   segment.hpp(\*)
compare_lines.cpp  io_png.h    mulsd.cpp(\*)  test_pyramid.cpp
(\*) Reviewed in IPOL.

[^1]: Rafael Grompone von Gioi, Jérémie Jakubowicz, Jean-Michel Morel, and Gregory Randall, LSD: a Line Segment Detector, Image Processing On Line, 2 (2012), pp. 35–55. <http://doi.org/10.5201/ipol.2012.gjmr-lsd>

[^2]: Yohann Salaün, Renaud Marlet, and Pascal Monasse, Multiscale line segment detector for robust and accurate SfM, Procedings of ICPR 2016. <https://doi.org/10.1109/ICPR.2016.7899930>

[^3]: Yohann Salaün, Renaud Marlet, and Pascal Monasse, The multiscale line segment detector, Proceedings of RRPR 2016. <https://doi.org/10.1007/978-3-319-56414-2_12>
