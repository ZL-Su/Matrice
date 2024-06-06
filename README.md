
# Geometric Optical Sensing and Visual Intelligence (R) 
# Matrice Library 
<img src="https://img.shields.io/github/v/release/ZL-Su/Matrice?include_prereleases&label=version"/>

A performancy primitives library for 2/3-dimensional visual perception and understanding and photomechanics with modern C++ programming language. The library, referred to as Matrice, covers the fundamental differentiable algebraic containers (including Scalar, Vector, Matrix and Tensor) and the associated algorithms (e.g., linear algebra operations), expression system, geometric algorithms, nonlinear optimization, and network building blocks.

## Features
* Unified memory and plain type management and runtime visualization.
* Automatic inference of optimal type and memory pattern for numeric containers according to the context.
* Powerful expression system with lazy evaluation, allowing for building efficient code with static or dynamic graph.
* Elegant and friendly algorithm implementations for linear systems, geomeric vision, feature recognition.

## Visualizer
Matrice supports runtime visualization for plain types, including dgelom::Matrix_<>, dgelom::Vec_<> and dgelom::tensor<>. At present, the Matrice visualizer is a extended script of the extension, Image Watch, in Visual Studio 2017/2019.

How-to:
* Install Image Watch from the online visual studio marketplace;
* Copy the visualizer script from "Matrice/Visualizers" to the Visualizer folder of Visual Studio IDE on local computer;
* Enjoy it.

## System Requirements
Matrice supports Intel x64 architecture and compatible architectures.
The library is optimized for the systems based on
* Intel(R) Core(TM) processor - 6th, 7th, 8th, 10th, and 11th generation
* Intel(R) Xeon(R) processor E5 v3 family (formerly Haswell)

and compatible processors.

The software dependencies are:
* C++ compiler with C++17 standart support (at least)
* Optional dependencies:
  * Intel OpenMP
  * CUDA 9.0 or later
## Reference
### *BibTex*
@misc{Su:2021,
  author = {Zhilong Su},
  year = {2021},
  url = {https://github.com/ZL-Su/Matrice},
  urldate = {February 02, 2020},
  title = {Matrice: A library of performance primitives in 3D vision and photomechanics}
}
