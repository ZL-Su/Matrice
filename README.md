<img src="../master/ivm.svg" width="800">

# IVM(R) Matrice Library 
<img src="../master/version.svg" width="80"/>

A libary for numerical computation, includes basic linear algebra (vector, matrix, linear system) operations, auto-differentiation, and non-linear optimization, etc. The library contains two parts, the first part is a Fortran kernel library which supplies the basic linear arithmatic operations; the second part is a modern C++ library which aims to implement efficient and elegant numeric computation.

## Features
* Unified memory and basic type management.
* Automatic inference of optimal type according to the context.
* Expression system enhanced lazy evaluation.

## Visualizer
Matrice supports runtime visualization for plain types, including Matrix_<>, Vec_<> and tensor<>. At present, the Matrice visualizer is a script for the extension, Image Watch, in Visual Studio 2017/2019.

How-to:
* Install Image Watch from the online visual studio marketplace;
* Copy the visualizer script from "Matrice/Visualizers" to where the Visualizer folder of Visual Studio IDE on local computer;
* Enjoy it.

## System Requirements
Matrice supports Intel 64 architecture and compatible architectures.
The library is optimized for the systems based on
* Intel(R) Core(TM) processor - 5th, 6th, 7th, and 8th generation
* Intel(R) Xeon(R) processor E5 v3 family (formerly Haswell)

and compatible processors.

The software dependencies are:
* C++ compiler with C++11 standart support (at least)
* Optinal dependencies:
  * Intel OpenMP
  * CUDA 8.0 or later
## Reference
### *BibTex*
@misc{Su:2020,
  author = {Zhilong Su},
  year = {2020},
  url = {https://github.com/ZL-Su/Matrice},
  urldate = {February 02, 2020},
  title = {Matrice Library}
}
