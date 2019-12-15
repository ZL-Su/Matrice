/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include "_tag_defs.h"
#include "_type_traits.h"
#include "storage/forward.hpp"

#if MATRICE_SIMD_ARCH == MATRICE_SIMD_AVX
#define MATRICE_ALIGN_BYTES   0x0020
#else
#ifdef MATRICE_SIMD_ARCH == MATRICE_SIMD_AVX512
#define MATRICE_ALIGN_BYTES   0x0040
#else
#define MATRICE_ALIGN_BYTES   0x0010
#endif
#endif

DGE_MATRICE_BEGIN
#ifdef MATRICE_ALIGN_BYTES
#define MATRICE_ALIGNED(TYPE) alignas(MATRICE_ALIGN_BYTES)##TYPE
#endif

enum Location {
	UnSpecified = -1,
	OnStack = 0,
	OnHeap = 1,
#ifdef MATRICE_ENABLE_CUDA
	OnDevice = 2,
	OnGlobal = 3,
#endif
};
enum {
	COPY = 1001,
	MOVE = 1002,
	SHARED = 1000
};
enum {
	LINEAR = 8000,
#ifdef MATRICE_ENABLE_CUDA
	PITCHED = 8001,
	ARRTARR = 8002,
	FROMARR = 8003,
	TOARRAY = 8004,
#endif
};

using loctn_t = Location; using memloc_t = Location;

namespace privt {
template<typename _Vty, typename _Ity> _Vty* aligned_malloc(_Ity size);
template<typename _Vty> void aligned_free(_Vty* aligned_ptr) noexcept;
template<typename _Vty> bool is_aligned(_Vty* aligned_ptr) noexcept;
template<typename _Vty, typename _Ity> MATRICE_HOST_FINL
constexpr _Vty* fill_mem(const _Vty* src, _Vty* dst, _Ity size)
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

#ifdef MATRICE_ENABLE_CUDA
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
_Ty* device_malloc(size_t& w, size_t h);
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
_Ty* global_malloc(size_t size);
template<size_t _Kind, typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
void device_memcpy(_Ty* hptr, _Ty* dptr, size_t w, size_t h = 1, size_t p = 1);
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
void device_free(_Ty* dptr);
template<int _Opt = 0> void _Device_sync();
#endif

}
namespace internal {
struct _Memory {
	/**
	 *\brief Check if a given memory is aligned or not
	 *\param [_Ptr] the pointer to memory block to be checked
	 *\param [_Aligns] align bytes for checking
	 */
	template<typename _It, MATRICE_ENABLE_IF(is_pointer_v<_It>)>
	inline static bool is_aligned(const _It _Ptr, size_t _Aligns = MATRICE_ALIGN_BYTES) {
		return !(reinterpret_cast<size_t>(reinterpret_cast<void*>(_Ptr))%_Aligns);
	}

	/**
	 *\brief Memory release
	 *\param [_Ptr] memory pointer
	 */
	template<Location _Loc, typename _It, MATRICE_ENABLE_IF(is_pointer_v<_It>)>
	MATRICE_HOST_INL static void free(_It _Ptr) {
		try {
			if /*constexpr */(_Loc == Location::OnHeap)
				if (is_aligned(_Ptr)) privt::aligned_free(_Ptr);
				else { std::free(_Ptr); _Ptr = nullptr; }
#ifdef MATRICE_ENABLE_CUDA
			else if /*constexpr*/ (_Loc == Location::OnGlobal)
				privt::device_free(_Ptr);
			else if /*constexpr*/ (_Loc == Location::OnDevice)
				privt::device_free(_Ptr);
#endif // MATRICE_ENABLE_CUDA
			else return;
		}
		catch (std::exception e) { throw e; }
	}

template<typename _InIt, typename _OutIt> 
MATRICE_HOST_INL static _OutIt copy(_InIt _First, _InIt _Last, _OutIt _Dest) {
	
	if (std::is_trivially_assignable_v<_OutIt, _InIt>) {
		const char* const _First_ch = const_cast<const char*>(reinterpret_cast<const volatile char*>(_First));
		const char* const _Last_ch = const_cast<const char*>(reinterpret_cast<const volatile char *>(_Last));
		char* const _Dest_ch = const_cast<char*>(reinterpret_cast<volatile char*>(_Dest));
		const auto _Count = static_cast<size_t>(_Last_ch - _First_ch);

		::memmove(_Dest_ch, _First_ch, _Count);

		return (reinterpret_cast<_OutIt>(_Dest_ch + _Count));
	}
	else {
		for (; _First != _Last; ++_Dest, (void)++_First)
		{
			*_Dest = *_First;
		}

		return (_Dest);
	}
}

};
}

namespace internal {
template<typename _Ty, class _Tag>
MATRICE_GLOBAL_INL decltype(auto) malloc(size_t rows, size_t& cols, _Tag);

template<typename _Ty, class _Tag>
MATRICE_GLOBAL_INL void free(_Ty* data, _Tag);

template<typename _Ty>
MATRICE_HOST_INL decltype(auto) aligned_malloc(size_t size);

template<typename _Ty>
MATRICE_HOST_INL void aligned_free(_Ty* aligned_ptr) noexcept;

template<typename _Ty>
MATRICE_HOST_INL bool is_aligned(_Ty* aligned_ptr) noexcept;

template<typename _InIt, typename _OutIt, class _InTag, class _OutTag>
MATRICE_GLOBAL_INL _OutIt copy(_InIt _First, _InIt _Last, _OutIt _Dest, _InTag, _OutTag);

template<class _Tag>
MATRICE_GLOBAL_INL decltype(auto) make_alloc_deleter(_Tag);

namespace impl {
template<typename _Ty>
MATRICE_HOST_INL decltype(auto) _Malloc(size_t rows, size_t cols, heap_alloc_tag);

template<typename _Ty>
MATRICE_HOST_INL void _Free(_Ty* data, heap_alloc_tag);

MATRICE_HOST_INL decltype(auto) _Deleter(stack_alloc_tag) noexcept;
MATRICE_HOST_INL decltype(auto) _Deleter(heap_alloc_tag) noexcept;

#ifdef MATRICE_ENABLE_CUDA
template<typename _Ty>
MATRICE_GLOBAL _Ty* _Malloc(size_t rows, size_t& cols, device_alloc_tag);

template<typename _Ty>
MATRICE_GLOBAL _Ty* _Malloc(size_t rows, size_t cols, global_alloc_tag);

template<typename _Ty>
MATRICE_HOST_INL void _Free(_Ty* data, device_alloc_tag);

template<typename _Ty>
MATRICE_HOST_INL void _Free(_Ty* data, global_alloc_tag);

MATRICE_HOST_INL decltype(auto) _Deleter(device_alloc_tag);
MATRICE_HOST_INL decltype(auto) _Deleter(global_alloc_tag);
#endif

template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt _Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, stack_alloc_tag, stack_alloc_tag);

template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt _Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, stack_alloc_tag, heap_alloc_tag);

template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt _Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, heap_alloc_tag, stack_alloc_tag);

template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt _Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, heap_alloc_tag, heap_alloc_tag);
}
}
DGE_MATRICE_END