/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
*********************************************************************/
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

#if _HAS_CXX20
#define __cpp_concepts 202007L
#endif

#define MATRICE_SIMD_DISABLE 0
#define MATRICE_SIMD_SSE    2 //*\SIMD-128
#define MATRICE_SIMD_AVX    3 //*\SIMD-256
#define MATRICE_SIMD_AVX512 4 //*\SIMD-512

#define MATRICE_USE_NAT     0 //*\use native kernel code
#define MATRICE_USE_MKL     1 //*\use intel mkl lib

#ifdef MATRICE_USE_OPENCV
#ifndef __use_ocv_as_view__
#define __use_ocv_as_view__
#endif
#endif // __enable_ocv__

#ifdef MATRICE_ENABLE_CUDA
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

#define MATRICE_NAMESPACE_BEGIN_ namespace dgelom {
#define _MATRICE_NAMESPACE_END }

#define DGE_MATRICE_BEGIN namespace dgelom {
#define DGE_MATRICE_END }

#define MATRICE_EXPR_BEGIN DGE_MATRICE_BEGIN \
inline namespace exprs {
#define MATRICE_EXPR_END } DGE_MATRICE_END

#define MATRICE_DEVICE_BEGIN DGE_MATRICE_BEGIN \
namespace device {
#define MATRICE_DEVICE_END  } DGE_MATRICE_END

#define MATRICE_PRIVATE_BEGIN DGE_MATRICE_BEGIN \
namespace privt {
#define MATRICE_PRIVATE_END  } DGE_MATRICE_END

#define MATRICE_ALGS_BEGIN DGE_MATRICE_BEGIN \
namespace algs {
#define MATRICE_ALGS_END  } DGE_MATRICE_END

#define MATRICE_ARCH_BEGIN DGE_MATRICE_BEGIN \
namespace simd {
#define MATRICE_ARCH_END   } DGE_MATRICE_END

#define MATRICE_ALG_BEGIN(ALG) DGE_MATRICE_BEGIN \
namespace ALG {
#define MATRICE_ALG_END(ALG) } DGE_MATRICE_END

#define _CONDITIONS_BEGIN namespace conds {
#define _CONDITIONS_END }
#define _DETAIL_BEGIN namespace detail {
#define _DETAIL_END }
#define _INTERNAL_BEGIN namespace internal {
#define _INTERNAL_END }

#ifndef MATRICE_CONSTEXPR_IF
#if defined (_HAS_CXX17)
#define MATRICE_CONSTEXPR_IF(_COND) \
 if constexpr (_COND)
#else
#define MATRICE_CONSTEXPR_IF(_COND) \
 if           (_COND)
#endif
#endif
#ifndef MATRICE_ENABLE_IF
#define MATRICE_ENABLE_IF(_COND) \
typename = std::enable_if_t<_COND>
#endif

#ifndef MATRICE_DEBUG
#ifdef _DEBUG
#define MATRICE_DEBUG
#endif
#endif

#ifdef _MSC_VER
#ifndef _ENABLE_EXTENDED_ALIGNED_STORAGE
#define _ENABLE_EXTENDED_ALIGNED_STORAGE
#endif
#endif

#ifndef MATRICE_USE_STD
#define MATRICE_USE_STD(NAME) \
using std::NAME;
#endif