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
#define MATRICE_SHARED_STORAGE 0
#endif //enable shared memory allocator

#ifndef __enable_cuda__ 
#define __enable_cuda__
#define MATRICE_ENABLE_CUDA __enable_cuda__
#endif // !__enable_cuda__ defined for CUDA support

#ifdef _OPENMP
#define __enable_omp__
#define MATRICE_ENABLE_OMP
#endif // _OPENMP


#ifndef __AVX__
#define __AVX__
#endif // !__AVX__
#ifndef __SSE__
#define __SSE__
#endif // !__SSE__

#define MATRICE_SIMD_SSE    2 //*\SIMD-128
#define MATRICE_SIMD_AVX    3 //*\SIMD-256
#define MATRICE_SIMD_AVX512 4 //*\SIMD-512

#ifndef MATRICE_SIMD_ARCH
#define MATRICE_SIMD_ARCH MATRICE_SIMD_AVX
#endif

#define MATRICE_USE_NAT     0 //*\use native kernel code
#define MATRICE_USE_FKL     1 //*\use fortran kernel lib
#define MATRICE_USE_MKL     2 //*\use intel mkl lib
#ifndef MATRICE_MATH_KERNEL
#define MATRICE_MATH_KERNEL MATRICE_USE_MKL
#endif // !MATRICE_MATH_KERNEL


#ifdef MATRICE_USE_OPENCV
#ifndef __use_ocv_as_view__
#define __use_ocv_as_view__
#endif
#endif // __enable_ocv__

#if (defined MATRICE_ENABLE_CUDA && !defined MATRICE_DISABLE_CUDA)
#define __disable_simd__
#define MATRICE_DISABLE_SIMD __disable_simd__
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

#ifndef MATRICE_NONHERITABLE
#define MATRICE_NONHERITABLE final
#endif

#ifndef __
#define __ 0      //run-time deduced on host
#endif

#ifndef _RTDD
#define _RTDD -1  //run-time deduced on device
#endif

#ifndef MATRICE_INT
#define MATRICE_INT                                std::ptrdiff_t
typedef MATRICE_INT                                     int_t;
#else
typedef int                                             int_t;
#endif

#define MATRICE_NAMESPACE_BEGIN_ namespace dgelom {
#define _MATRICE_NAMESPACE_END                   }
#define DGE_MATRICE_BEGIN namespace dgelom {
#define DGE_MATRICE_END }
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
#define MATRICE_ARCH_BEGIN MATRICE_NAMESPACE_BEGIN_ namespace simd {
#define MATRICE_ARCH_END   } _MATRICE_NAMESPACE_END
#define _TYPES_BEGIN namespace types {
#define _TYPES_END }
#define _CONDITIONS_BEGIN namespace conds {
#define _CONDITIONS_END }
#define _DETAIL_BEGIN namespace detail {
#define _DETAIL_END }
#define _INTERNAL_BEGIN namespace internal {
#define _INTERNAL_END }

#ifndef MATRICE_CONSTEXPR_IF
#if defined (_HAS_CXX17)
#define MATRICE_CONSTEXPR_IF(_COND) if constexpr (_COND)
#else
#define MATRICE_CONSTEXPR_IF(_COND) if           (_COND)
#endif
#endif
#ifndef MATRICE_ENABLE_IF
#define MATRICE_ENABLE_IF(_COND) typename = std::enable_if_t<_COND>
#endif