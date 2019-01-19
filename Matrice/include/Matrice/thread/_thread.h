/*******************************************************************************
* Copyright 2018-2019 IVM (R) DGELOM
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*******************************************************************************/
#pragma once

#include <type_traits>

#define MATRICE_THR_SEQ 0
#define MATRICE_THR_OMP 1
#define MATRICE_THR_TBB 2

/* Ideally this condition below should never happen (if the library is built
 * using regular cmake). For the 3rd-party projects that build the library
 * from the sources on their own try to guess the right threading... */
#if !defined(MATRICE_THR)
#   if defined(_OPENMP)
#       define MATRICE_THR MATRICE_THR_OMP
#   else
#       define MATRICE_THR MATRICE_THR_SEQ
#   endif
#endif

#define CHAIn2(a,b) a b
#define CHAIN2(a,b) CHAIn2(a,b)

#define CONCAt2(a,b) a ## b
#define CONCAT2(a,b) CONCAt2(a,b)

#define STRINGIFy(s) #s
#define STRINGIFY(s) STRINGIFy(s)

#ifdef _MSC_VER
#   define PRAGMA_MACRo(x) __pragma(x)
#   define PRAGMA_MACRO(x) PRAGMA_MACRo(x)
#else
#   define PRAGMA_MACRo(x) _Pragma(#x)
#   define PRAGMA_MACRO(x) PRAGMA_MACRo(x)
#endif

#if MATRICE_THR == MATRICE_THR_SEQ
#define MATRICE_THR_SYNC 1
inline int _Get_max_threads() { return 1; }
inline int _Get_num_threads() { return 1; }
inline int _Get_thread_num() { return 0; }
inline int _In_parallel() { return 0; }
inline void _Host_thr_barrier() {}

#define PRAGMA_OMP(...)

#elif MATRICE_THR == MATRICE_THR_OMP
#include <omp.h>
#define MATRICE_THR_SYNC 1

inline int _Get_max_threads() { return omp_get_max_threads(); }
inline int _Get_num_threads() { return omp_get_num_threads(); }
inline int _Get_thread_num() { return omp_get_thread_num(); }
inline int _In_parallel() { return omp_in_parallel(); }
inline void _Host_thr_barrier() {
#   pragma omp barrier
}

#define PRAGMA_OMP(...) PRAGMA_MACRO(CHAIN2(omp, __VA_ARGS__))

#elif MATRICE_THR == MATRICE_THR_TBB
#include "tbb/task_arena.h"
#include "tbb/parallel_for.h"
#define MATRICE_THR_SYNC 0

inline int _Get_max_threads()
{ return tbb::this_task_arena::max_concurrency(); }
inline int _Get_num_threads() { return _Get_max_threads(); }
inline int _Get_thread_num()
{ return tbb::this_task_arena::current_thread_index(); }
inline int _In_parallel() { return 0; }
inline void _Host_thr_barrier() { assert(!"no barrier in TBB"); }

#define PRAGMA_OMP(...)

#endif

/* MSVC still supports omp 2.0 only */
#if defined(_MSC_VER) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#   define collapse(x)
#   define PRAGMA_OMP_SIMD(...)
#else
#   define PRAGMA_OMP_SIMD(...) PRAGMA_MACRO(CHAIN2(omp, simd __VA_ARGS__))
#endif // defined(_MSC_VER) && !defined(__INTEL_COMPILER)

namespace dgelom {
namespace impl {

inline bool _Host_thr_syncable() { return MATRICE_THR_SYNC == 1; }

namespace utils {
	template <typename T, typename U>
	inline std::remove_reference_t<T> div_up(const T a, const U b) {
		assert(b);
		return (a + b - 1) / b;
	}

	template <typename T, typename U>
	inline std::remove_reference_t<T> rnd_up(const T a, const U b) {
		return div_up(a, b) * b;
	}

	template <typename T, typename U>
	inline std::remove_reference_t<T> rnd_dn(const T a, const U b) {
		return (a / b) * b;
	}
}

template<size_t _Nt, typename _Fn> inline void parallel(_Fn&& fn);
template<typename... _Args> inline void for_nd(const int ithr, const int _Nt, const _Args&... _args);
template<typename... _Args> inline void parallel_nd(_Args&&... args);
} // namespace impl
using impl::parallel;
using impl::parallel_nd;
} // namespace dgelom

#include "inline\_thread_parallel_nd.inl"

// vim: et ts=4 sw=4 cindent cino^=l0,\:0,N-s
