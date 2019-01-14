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

#include <xutility>
#include "..\util\_macros.h"
#include "_tag_defs.h"

#ifdef __AVX__
#define MATRICE_ALIGN_BYTES   0x0020
#else
#ifdef __AVX512
#define MATRICE_ALIGN_BYTES   0x0040
#else
#define MATRICE_ALIGN_BYTES   0x0010
#endif
#endif // __AVX__

DGE_MATRICE_BEGIN
#ifdef MATRICE_ALIGN_BYTES
#define MATRICE_ALIGNED(TYPE) alignas(MATRICE_ALIGN_BYTES)##TYPE
#endif

enum Location { UnSpecified = -1, OnStack = 0, OnHeap = 1, OnDevice = 2, OnGlobal = 3 };
using loctn_t = Location; using memloc_t = Location;

enum { COPY = 1001, MOVE = 1002, SHARED = 1000 };
enum { LINEAR = 8000, PITCHED = 8001, ARRTARR = 8002, FROMARR = 8003, TOARRAY = 8004, };

namespace privt {
template<typename ValueType, typename IntegerType> ValueType* aligned_malloc(IntegerType size);
template<typename ValueType> void aligned_free(ValueType* aligned_ptr) noexcept;
template<typename ValueType> bool is_aligned(ValueType* aligned_ptr) noexcept;
template<typename ValueType, typename Integer> MATRICE_HOST_FINL
ValueType* fill_mem(const ValueType* src, ValueType* dst, Integer size)
{
#define _FILOP(_Idx) dst[_Idx] = src[_Idx]
#define _RET return (dst);

	if (size == 1) { _FILOP(0); _RET; }
	if (size == 2) { _FILOP(0), _FILOP(1); _RET; }
	if (size == 3) { _FILOP(0), _FILOP(1), _FILOP(2); _RET; }
	if (size == 4) { _FILOP(0), _FILOP(1), _FILOP(2), _FILOP(3); _RET; }
	if (size == 5) { _FILOP(0), _FILOP(1), _FILOP(2), _FILOP(3), _FILOP(4); _RET; }
	if (size == 6) { _FILOP(0), _FILOP(1), _FILOP(2), _FILOP(3), _FILOP(4), _FILOP(5); _RET; }
	if (size == 7) { _FILOP(0), _FILOP(1), _FILOP(2), _FILOP(3), _FILOP(4), _FILOP(5), _FILOP(6); _RET; }
	if (size == 8) { _FILOP(0), _FILOP(1), _FILOP(2), _FILOP(3), _FILOP(4), _FILOP(5), _FILOP(6), _FILOP(7); _RET; }
	if (size == 9) { _FILOP(0), _FILOP(1), _FILOP(2), _FILOP(3), _FILOP(4), _FILOP(5), _FILOP(6), _FILOP(7), _FILOP(8); _RET; }

	if (size > 9) { std::copy(src, src + size, dst); _RET; }

#undef _FILOP
#undef _RET
}
}
namespace internal {
struct _Memory {

	/**
	 *\brief Check if a given memory is aligned or not
	 *\param [_Ptr] the pointer to memory block to be checked
	 *\param [_Aligns] align bytes for checking
	 */
	template<typename _It>
	static bool is_aligned(const _It _Ptr, size_t _Aligns = MATRICE_ALIGN_BYTES) {
		return !(reinterpret_cast<size_t>(reinterpret_cast<void*>(_Ptr)) % _Aligns);
	}

	/**
	 *\brief Memory release
	 *\param [_Ptr] memory pointer
	 */
	template<typename _It, 
		typename = std::enable_if_t<std::is_pointer_v<_It>>>
	static void free(_It _Ptr) {
		try {
			if (is_aligned(_Ptr)) privt::aligned_free(_Ptr);
			else std::free(_Ptr);
		}
		catch (std::exception e) {}
	}
};
}
DGE_MATRICE_END