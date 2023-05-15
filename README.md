# LR B-splines for matlab


## Introduction

LR B-splines is a technology which enables users to perform local refinement on B-spline surfaces or volumes, which traditionally have been limited to tensor products. This has a wide variety of applications within computer-aided design (CAD) and computer-aided engineering (CAE). While this library was written with the latter in mind, it is also possible to take use of it in a design environment.

## About the code

This is a matlab wrapper over the core c++ library. You will have to compile the c++ library first, followed by linking the `.mex` files provided here to this library. This will in turn allow give you the full power of a fast c++ library with the convenient matlab syntax on top.

## Compiling on Ubuntu

Assuming you have installed [LR B-splines](https://github.com/VikingScientist/LRsplines) and have a version of [Matlab](https://se.mathworks.com/products/matlab.html) installed, then the library is compiled and linked by simply typing
```
cmake .
make
```
in the root folder.

## Compiling in Matlab
Copy all `.cpp` and `.h` files from [LR B-splines](https://github.com/VikingScientist/LRsplines). It is sufficient to copy over the entire `src/`- and `include/`-folders. After this add the Eigen3 library in the folder `eigen3` in this project root folder. Eigen can be downloaded from the [eigen3 website](https://eigen.tuxfamily.org/index.php?title=Main_Page). The library can then be compiled by the command
```
mkdir lib
mex src/*.cpp -Iinclude -output lib/lrsplinesurface_interface
copyfile src/matlab/* lib
```
To start using the library, add this to the path
```
addpath('lib')
```
And test it all by running `lrgui`.
