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

#ifdef __AVX__
#define MATRICE_ALIGN_BYTES   0x0020
#else
#ifdef __AVX512
#define MATRICE_ALIGN_BYTES   0x0040
#else
#define MATRICE_ALIGN_BYTES   0x0010
#endif
#endif // __AVX__

namespace dgelom {
#ifdef MATRICE_ALIGN_BYTES
#define MATRICE_ALIGNED(type) alignas(MATRICE_ALIGN_BYTES)##type
#endif
	typedef enum Location
	{
		UnSpecified = -1, OnStack = 0, OnHeap = 1, OnDevice = 2, OnGlobal = 3
	} loctn_t;
	enum { COPY = 1001, MOVE = 1002, SHARED = 1000 };
	enum { LINEAR = 8000, PITCHED = 8001, ARRTARR = 8002, FROMARR = 8003, TOARRAY = 8004, };
namespace privt {
template<typename ValueType, typename IntegerType> ValueType* aligned_malloc(IntegerType size);
template<typename ValueType> void aligned_free(ValueType* aligned_ptr) noexcept;
template<typename ValueType> bool is_aligned(ValueType* aligned_ptr) noexcept;
template<typename ValueType, typename Integer> __forceinline 
ValueType* fill_mem(const ValueType* src, ValueType* dst, Integer size)
{
	if (size == 1) 
		dst[0] = src[0];
	if (size == 2)
		dst[0] = src[0], dst[1] = src[1];
	if (size == 3)
		dst[0] = src[0], dst[1] = src[1], dst[2] = src[2];
	if (size == 4)
		dst[0] = src[0], dst[1] = src[1], dst[2] = src[2], dst[3] = src[3];
	if (size == 5)
		dst[0] = src[0], dst[1] = src[1], dst[2] = src[2], dst[3] = src[3], dst[4] = src[4];
	if (size == 6)
		dst[0] = src[0], dst[1] = src[1], dst[2] = src[2], dst[3] = src[3], dst[4] = src[4], dst[5] = src[5];
	if (size == 7)
		dst[0] = src[0], dst[1] = src[1], dst[2] = src[2], dst[3] = src[3], dst[4] = src[4], dst[5] = src[5], dst[6] = src[6];
	if (size == 8)
		dst[0] = src[0], dst[1] = src[1], dst[2] = src[2], dst[3] = src[3], dst[4] = src[4], dst[5] = src[5], dst[6] = src[6], dst[7] = src[7];
	if (size == 9)
		dst[0] = src[0], dst[1] = src[1], dst[2] = src[2], dst[3] = src[3], dst[4] = src[4], dst[5] = src[5], dst[6] = src[6], dst[7] = src[7], dst[8] = src[8];
	if (size > 9) std::copy(src, src + size, dst);
	return (dst);
}



}
}