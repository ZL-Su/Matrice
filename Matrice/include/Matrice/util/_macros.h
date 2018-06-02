/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once

#ifndef __CXX11__
#define __CXX11__
#endif // !__CXX11__ defined for C++11 standard support

#ifndef __CXX17__
#define __CXX17__
#endif // !__CXX17__ defined for C++17 standard support

#if (defined __CXX11__ || defined __CXX17__)
#define __CXX11_SHARED__
#endif //enable shared memory allocator

#ifndef __enable_cuda__
#define __enable_cuda__
#endif // !__enable_cuda__ defined for CUDA support

#ifdef _OPENMP
#define __enable_omp__
#endif // _OPENMP


#ifndef __AVX__
#define __AVX__
#endif // !__AVX__
#ifndef __SSE__
#define __SSE__
#endif // !__SSE__

#ifdef __enable_ocv__
#ifndef __use_ocv_as_view__
#define __use_ocv_as_view__
#endif
#endif // __enable_ocv__

#if (defined __enable_cuda__ && !defined __disable_cuda__)
#define __disable_simd__
#include <host_defines.h>
#define MATRICE_HOST_ONLY __host__
#define MATRICE_DEVICE_ONLY __device__
#define MATRICE_GLOBAL __host__ __device__
#define MATRICE_HOST_INL __inline__ __host__
#define MATRICE_DEVICE_INL __inline__ __device__
#define MATRICE_GLOBAL_INL __inline__ __host__ __device__
#define MATRICE_HOST_FINL __forceinline__ __host__
#define MATRICE_DEVICE_FINL __forceinline__ __device__
#define MATRICE_GLOBAL_FINL __forceinline__ __host__ __device__
#else
#define MATRICE_HOST_ONLY
#define MATRICE_DEVICE_ONLY 
#define MATRICE_GLOBAL
#define MATRICE_HOST_INL inline
#define MATRICE_DEVICE_INL inline
#define MATRICE_GLOBAL_INL __inline
#define MATRICE_HOST_FINL __forceinline
#define MATRICE_DEVICE_FINL __forceinline
#define MATRICE_GLOBAL_FINL __forceinline
#endif

#ifndef __
#define __ 0
#endif

#ifndef _1
#define _1 -1
#endif

#define MATRICE_NAMESPACE_BEGIN_ namespace dgelom {
#define _MATRICE_NAMESPACE_END                   }
#define MATRICE_NAMESPACE_BEGIN_TYPES MATRICE_NAMESPACE_BEGIN_ namespace types {
#define MATRICE_NAMESPACE_END_TYPES  } _MATRICE_NAMESPACE_END
#define MATRICE_NAMESPACE_EXPR_BEGIN MATRICE_NAMESPACE_BEGIN_ namespace exprs {
#define MATRICE_NAMESPACE_EXPR_END MATRICE_NAMESPACE_END_TYPES
#define MATRICE_DEVICE_BEGIN MATRICE_NAMESPACE_BEGIN_ namespace device {
#define MATRICE_DEVICE_END  } _MATRICE_NAMESPACE_END
#define MATRICE_PRIVATE_BEGIN MATRICE_NAMESPACE_BEGIN_ namespace privt {
#define MATRICE_PRIVATE_END  } _MATRICE_NAMESPACE_END
#define MATRICE_ALGS_BEGIN MATRICE_NAMESPACE_BEGIN_ namespace algs {
#define MATRICE_ALGS_END  } _MATRICE_NAMESPACE_END