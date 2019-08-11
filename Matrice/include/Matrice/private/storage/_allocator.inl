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

#define MATRICE_MEMCPY_ADAPTER_1 \
typename _Al::category(), category()
#define MATRICE_MEMCPY_ADAPTER_2 \
category(), typename _Al::category()

DGE_MATRICE_BEGIN
_DETAIL_BEGIN

template<typename _Altrs> template<typename _Al>
MATRICE_GLOBAL_INL _Dense_allocator_base<_Altrs>& _Dense_allocator_base<_Altrs>::_Alloc_copy(const _Al& al) noexcept {
	const auto _First = al.data();
	const auto _Last = _First + min(size(), al.size());
	internal::copy(_First, _Last, data(), MATRICE_MEMCPY_ADAPTER_1);
	return (*this);
}

template<typename _Altrs> template<typename _Al>
MATRICE_GLOBAL_INL _Dense_allocator_base<_Altrs>& _Dense_allocator_base<_Altrs>::_Alloc_move(_Al&& al) noexcept {
	m_data = al.data();
	m_rows = al.rows();
	m_cols = al.cols();
	al.data() = nullptr;
	return (*this);
}

MATRICE_ALLOCATOR(HOST,, 0)::~_Allocator() {
	internal::aligned_free(this->m_data);
}

MATRICE_ALLOCATOR(HOST, decltype(auto), 0)::_Alloc() noexcept {
	this->data() = internal::aligned_malloc<value_type>(this->size());
	return (*this);
}
_DETAIL_END

namespace internal {
template<typename _Ty>
MATRICE_HOST_INL decltype(auto) aligned_malloc(size_t size) {
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
MATRICE_HOST_INL void aligned_free(_Ty* aligned_ptr) noexcept {
	if (aligned_ptr) {
		std::free(*(reinterpret_cast<void**>(reinterpret_cast<void*>(aligned_ptr)) - 1));
		aligned_ptr = nullptr;
	}
}

template<typename _Ty>
MATRICE_HOST_INL bool is_aligned(_Ty* aligned_ptr) noexcept {
	return !(reinterpret_cast<size_t>(reinterpret_cast<void*>(aligned_ptr)) % MATRICE_ALIGN_BYTES);
}

template<typename _InIt, typename _OutIt, class _InCat, class _OutCat>
MATRICE_GLOBAL_INL _OutIt copy(_InIt _First, _InIt _Last, _OutIt _Dest, _InCat, _OutCat) {
	return impl::_Memcpy(_First, _Last, _Dest, _InCat(), _OutCat());
}

namespace impl {
template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt _Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, stack_alloc_tag, stack_alloc_tag) {
	for (; _First != _Last; ++_Dest, (void)++_First) {
		*_Dest = *_First;
	}
	return (_Dest);
}

template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt _Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, stack_alloc_tag, heap_alloc_tag) {
	for (; _First != _Last; ++_Dest, (void)++_First) {
		*_Dest = *_First;
	}
	return (_Dest);
}

template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt _Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, heap_alloc_tag, stack_alloc_tag) {
	for (; _First != _Last; ++_Dest, (void)++_First) {
		*_Dest = *_First;
	}
	return (_Dest);
}

template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt _Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, heap_alloc_tag, heap_alloc_tag) {
	if (std::is_trivially_assignable_v<_OutIt, _InIt> || is_same_v<_InIt, _OutIt>) {
		const char* const _First_ch = const_cast<const char*>(reinterpret_cast<const volatile char*>(_First));
		const char* const _Last_ch = const_cast<const char*>(reinterpret_cast<const volatile char *>(_Last));
		char* const _Dest_ch = const_cast<char*>(reinterpret_cast<volatile char*>(_Dest));
		const auto _Count = static_cast<size_t>(_Last_ch - _First_ch);

		::memmove(_Dest_ch, _First_ch, _Count);

		return (reinterpret_cast<_OutIt>(_Dest_ch + _Count));
	}
	else {
		for (; _First != _Last; ++_Dest, (void)++_First) {
			*_Dest = *_First;
		}
		return (_Dest);
	}
}
}
}
DGE_MATRICE_END

#undef MATRICE_ALLOCATOR
#undef MATRICE_MEMCPY_ADAPTER_1
#undef MATRICE_MEMCPY_ADAPTER_2