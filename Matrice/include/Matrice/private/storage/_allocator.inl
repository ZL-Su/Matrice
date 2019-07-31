/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for
more detail.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once

#include "../_storage.hpp"

#define MATRICE_ALLOCATOR(LOC, RET, SV) \
template<typename _Ty, typename _Layout> \
 \
MATRICE_##LOC##_INL RET \
_Allocator<_Ty, SV, SV, \
allocator_traits_v<SV,SV>, \
_Layout> 

DGE_MATRICE_BEGIN
namespace internal {
	template<typename _Ty>
	decltype(auto) aligned_malloc(size_t size) {
		using value_type = _Ty;
		try {
			auto raw_ptr = std::malloc(size * sizeof(value_type) + MATRICE_ALIGN_BYTES);
			auto space = reinterpret_cast<size_t>(raw_ptr);
			space = space & ~(size_t(MATRICE_ALIGN_BYTES - 1));
			auto aligned_ptr = reinterpret_cast<void*>(space + MATRICE_ALIGN_BYTES);
			*(reinterpret_cast<void**>(aligned_ptr) - 1) = raw_ptr;

			return (reinterpret_cast<value_type*>(aligned_ptr));
		}
		catch (std::bad_alloc) {
			std::exception("Bad memory allocation.");
		};
	}
	template<typename _Ty>
	void aligned_free(_Ty* aligned_ptr) noexcept {
		if (aligned_ptr) {
			std::free(*(reinterpret_cast<void**>(reinterpret_cast<void*>(aligned_ptr)) - 1));
			aligned_ptr = nullptr;
		}
	}
	template<typename _Ty>
	bool is_aligned(_Ty* aligned_ptr) noexcept {
		return !(reinterpret_cast<size_t>(reinterpret_cast<void*>(aligned_ptr)) % MATRICE_ALIGN_BYTES);
	}
}

_DETAIL_BEGIN

MATRICE_ALLOCATOR(HOST, decltype(auto), 0)::_Alloc() noexcept {
	return internal::aligned_malloc<value_type>(this->size());
}

MATRICE_ALLOCATOR(HOST,, 0)::~_Allocator() {
	internal::aligned_free(this->m_data);
}

_DETAIL_END
DGE_MATRICE_END

#undef MATRICE_ALLOCATOR