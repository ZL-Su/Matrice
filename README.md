![Optinal Text](../master/ivm_logo_hori.png)

# IVM(R) Matrice Library

A libary for numerical computation, includes basic linear algebra (vector, matrix, linear system) operations, auto-differential, and non-linear optimization, etc. The library contains two parts, the first part is a Fortran kernel library which supplies the basic linear arithmatic operations; the second part is a modern C++ library which aims to implement efficient and elegant numeric computation.

## System Requirements
Matrice supports Intel 64 architecture and compatible architectures.
The library is optimized for the systems based on
* Intel(R) Core(TM) processor - 5th, 6th, 7th, and 8th generation
* Intel(R) Xeon(R) processor E5 v3 family (formerly Haswell)

and compatible processors.

The software dependencies are:
* C++ compiler with C++11 standart support
* Optinal dependencies:
  * Intel OpenMP
  * CUDA 8.0 or later
