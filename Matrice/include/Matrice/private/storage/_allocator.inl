/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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
#include "../_memory.h"

#define MATRICE_ALLOCTOR_SIG(M, N, OPT) \
_Allocator<_Ty, M, N, allocator_traits_v<OPT>, _Layout>

#define MATRICE_MEMCPY_ADAPTER_1 \
typename _Al::category(), category()
#define MATRICE_MEMCPY_ADAPTER_2 \
category(), typename _Al::category()

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Altrs>
MATRICE_GLOBAL_INL typename _Dense_allocator_base<_Altrs>::allocator&
_Dense_allocator_base<_Altrs>::operator=(const allocator& othr) noexcept {
	if (!(this->rows() == othr.rows() && this->cols()==othr.cols())) {
		alloc(othr.rows(), othr.cols());
	}
	if (this != &othr) {
		_Alloc_copy(othr.derived());
	}
	return (this->derived());
}

template<typename _Altrs>
MATRICE_GLOBAL_INL typename _Dense_allocator_base<_Altrs>::allocator&
_Dense_allocator_base<_Altrs>::operator=(allocator&& othr) noexcept {
	if (this != &othr) {
		_Alloc_move(othr.derived());
	}
	return (this->derived());
}

template<typename _Altrs>
MATRICE_GLOBAL_INL typename _Dense_allocator_base<_Altrs>::allocator&
_Dense_allocator_base<_Altrs>::operator=(const value_type val) noexcept {
	if (m_data != nullptr) {
		for (auto idx = 0; idx < this->size(); ++idx) {
			m_data[idx] = val;
		}
	}
	return (this->derived());
}

template<typename _Altrs>
MATRICE_GLOBAL_INL typename _Dense_allocator_base<_Altrs>::allocator&
_Dense_allocator_base<_Altrs>::operator=(const pointer data) noexcept {
	if (m_data != nullptr && data) {
		for (auto idx = 0; idx < this->size(); ++idx) {
			m_data[idx] = data[idx];
		}
	}
	return (this->derived());
}

template<typename _Altrs>
MATRICE_GLOBAL_INL void _Dense_allocator_base<_Altrs>::destroy() noexcept {
	if constexpr (is_same_v<category, heap_alloc_tag>) {
		internal::free(m_data, size(), category());
	}
#ifdef MATRICE_ENABLE_CUDA
	else if constexpr (is_any_of_v<category, device_alloc_tag, global_alloc_tag>) {
		internal::free(m_data, size(), category());
	}
#endif
}

template<typename _Altrs>
MATRICE_GLOBAL_INL decltype(auto) 
_Dense_allocator_base<_Altrs>::deleter() const noexcept {
	return internal::make_alloc_deleter(category());
}

template<typename _Altrs> template<typename _Al>
MATRICE_GLOBAL_INL typename _Dense_allocator_base<_Altrs>::allocator& 
_Dense_allocator_base<_Altrs>::_Alloc_copy(const _Al& al) noexcept {
	const auto _First = al.data();
	const auto _Last = _First + min(size(), al.size());
	internal::copy(_First, _Last, data(), MATRICE_MEMCPY_ADAPTER_1);
	return (this->derived());
}

template<typename _Altrs> template<typename _Al>
MATRICE_GLOBAL_INL typename _Dense_allocator_base<_Altrs>::allocator& 
_Dense_allocator_base<_Altrs>::_Alloc_move(_Al& al) noexcept {
	MATRICE_USE_STD(move);

	m_data = move(al.m_data);
	m_rows = move(al.m_rows);
	m_cols = move(al.m_cols);

	if constexpr (is_not_same_v<category, stack_alloc_tag>)
		al._Free_ownership();

	return (this->derived());
}

template<typename _Ty, typename _Layout>
MATRICE_HOST_INL decltype(auto)
MATRICE_ALLOCTOR_SIG(::dynamic, ::dynamic, ::dynamic)::_Alloc() noexcept {
	auto cols = _Mybase::cols();
	this->data() = internal::malloc<value_type>(_Mybase::rows(), cols, 
		typename _Mybase::category());
	return (*this);
}

//template<typename _Ty, typename _Layout>
//MATRICE_HOST_INL MATRICE_ALLOCTOR_SIG(::dynamic, ::dynamic, ::dynamic)::~_Allocator() {
//	if(this->size())
//	internal::free(this->data(), typename _Mybase::category());
//}

template<typename _Ty, diff_t _Rows, typename _Layout>
MATRICE_HOST_INL decltype(auto)
MATRICE_ALLOCTOR_SIG(_Rows, ::dynamic, ::dynamic)::_Alloc() noexcept {
	size_t cols = _Mybase::cols();
	this->data() = internal::malloc<value_type>(rows(), cols, 
		typename _Mybase::category());
	return (*this);
}

//template<typename _Ty, diff_t _Rows, typename _Layout>
//MATRICE_HOST_INL MATRICE_ALLOCTOR_SIG(_Rows, ::dynamic, ::dynamic)::~_Allocator() {
//	if (this->size())
//	internal::free(this->data(), typename _Mybase::category());
//}

template<typename _Ty, diff_t _Cols, typename _Layout>
MATRICE_HOST_INL decltype(auto)
MATRICE_ALLOCTOR_SIG(::dynamic, _Cols, ::dynamic)::_Alloc() noexcept {
	size_t cols = _Mybase::cols_at_compiletime;
	this->data() = internal::malloc<value_type>(_Mybase::rows(), cols, 
		typename _Mybase::category());
	return (*this);
}

//template<typename _Ty, diff_t _Cols, typename _Layout>
//MATRICE_HOST_INL MATRICE_ALLOCTOR_SIG(::dynamic, _Cols, ::dynamic)::~_Allocator() {
//	if (this->size())
//	internal::free(this->data(), typename _Mybase::category());
//}

#ifdef MATRICE_ENABLE_CUDA
template<typename _Ty, typename _Layout>
MATRICE_GLOBAL_INL decltype(auto)
MATRICE_ALLOCTOR_SIG(::device, ::device, ::device)::_Alloc() noexcept {
	m_pitch = this->cols();
	this->data() = internal::malloc<value_type>(this->rows(), m_pitch,
		device_alloc_tag());
	return (*this);
}

//template<typename _Ty, typename _Layout>
//MATRICE_GLOBAL_INL MATRICE_ALLOCTOR_SIG(::device, ::device, ::device)::~_Allocator() {
//	internal::free(this->data(), typename _Mybase::category());
//}

template<typename _Ty, typename _Layout>
MATRICE_GLOBAL_INL decltype(auto)
MATRICE_ALLOCTOR_SIG(::global, ::global, ::global)::_Alloc() noexcept {
	auto cols = this->cols();
	this->data() = internal::malloc<value_type>(this->rows(), cols, 
		global_alloc_tag());
	return (*this);
}

//template<typename _Ty, typename _Layout>
//MATRICE_GLOBAL_INL MATRICE_ALLOCTOR_SIG(::global, ::global, ::global)::~_Allocator() {
//	internal::free(this->data(), typename _Mybase::category());
//}
#endif
_DETAIL_END

namespace internal {
template<typename _Ty>
MATRICE_GLOBAL_INL constexpr const _Ty& _Maxval(const _Ty& left, const _Ty& right) noexcept {
	return left < right ? right : left;
}

template<typename _Ty, class _Tag>
MATRICE_GLOBAL_INL _Ty* malloc(size_t rows, size_t& cols, _Tag) {
	return impl::_Malloc<_Ty>(rows*cols, _Tag());
}

template<typename _Ty, class _Tag>
MATRICE_GLOBAL_INL void free(_Ty* data, _Tag) {
	return impl::_Free(data, _Tag());
}

template<typename _Ty, class _Tag>
void free(_Ty* data, size_t size, _Tag)
{
	return impl::_Free(data, size, _Tag());
}

template<typename _Ty>
MATRICE_HOST_INL _Ty* aligned_malloc(size_t size) {
	MATRICE_USE_STD(malloc);
	MATRICE_USE_STD(bad_alloc);
	try {
		auto raw_ptr = malloc(size * sizeof(_Ty) + MATRICE_ALIGN_BYTES);
		auto space = reinterpret_cast<size_t>(raw_ptr);
		space = space & ~(size_t(MATRICE_ALIGN_BYTES - 1));
		auto aligned_ptr = reinterpret_cast<void*>(space + MATRICE_ALIGN_BYTES);
		*(reinterpret_cast<void**>(aligned_ptr) - 1) = raw_ptr;

		return (reinterpret_cast<_Ty*>(aligned_ptr));
	}
	catch (bad_alloc) {
#ifdef MATRICE_DEBUG
		exception::error("Bad memory allocation in Func: aligned_malloc<_Ty>(size_t) in File: _allocator.inl.\n");
#else
		throw("Bad memory allocation in Func: aligned_malloc<_Ty>(size_t) in File: _allocator.inl.\n");
#endif
	};
	MATRICE_USE_STD(_Allocate);
	MATRICE_USE_STD(_New_alignof);
	MATRICE_USE_STD(_Get_size_of_n);
	return static_cast<_Ty*>(_Allocate<_New_alignof<_Ty>>(_Get_size_of_n<sizeof(_Ty)>(size)));
}

template<typename _Ty>
MATRICE_HOST_INL void aligned_free(_Ty* aligned_ptr) noexcept {
	MATRICE_USE_STD(free);
	if (aligned_ptr) {
		free(*(reinterpret_cast<void**>(reinterpret_cast<void*>(aligned_ptr)) - 1));
		aligned_ptr = nullptr;
	}
}

template<typename _Ty>
MATRICE_HOST_INL bool is_aligned(_Ty* aligned_ptr) noexcept {
	MATRICE_USE_STD(_New_alignof);
	return !(reinterpret_cast<size_t>(reinterpret_cast<void*>(aligned_ptr)) % 
		_Maxval<size_t>(_New_alignof<_Ty>, MATRICE_ALIGN_BYTES));
}

template<typename _InIt, typename _OutIt, class _InCat, class _OutCat>
MATRICE_GLOBAL_INL _OutIt copy(_InIt _First, _InIt _Last, _OutIt _Dest, _InCat, _OutCat) {
	return impl::_Memcpy(_First, _Last, _Dest, _InCat(), _OutCat());
}

template<typename _Tag>
MATRICE_GLOBAL_INL decltype(auto) make_alloc_deleter(_Tag) {
	return impl::_Deleter(_Tag());
}

namespace impl {
template<typename _Ty>
MATRICE_HOST_INL _Ty* _Malloc(size_t size, heap_alloc_tag) {
	MATRICE_USE_STD(_New_alignof);
	constexpr auto _Align = _Maxval<size_t>(_New_alignof<_Ty>, MATRICE_ALIGN_BYTES);
	return static_cast<_Ty*>(std::_Allocate<_Align>(size*sizeof(_Ty)));
}

template<typename _Ty>
MATRICE_HOST_INL void _Free(_Ty* data, heap_alloc_tag) {
	if (data != nullptr) {
		auto void_ptr = reinterpret_cast<void*>(data);
		if (void_ptr) {
			MATRICE_USE_STD(free);
			free(*(reinterpret_cast<void**>(void_ptr) - 1));
			data = nullptr;
		}
	}
}

template<typename _Ty>
MATRICE_HOST_INL void _Free(_Ty* data, size_t size, heap_alloc_tag) {
	constexpr size_t _Big_allocation_sentinel = 0xFAFAFAFAFAFAFAFAULL;
	MATRICE_USE_STD(_New_alignof);
	if (data && size < _Big_allocation_sentinel) {
		constexpr auto _Align = _Maxval<size_t>(_New_alignof<_Ty>, MATRICE_ALIGN_BYTES);
		std::_Deallocate<_Align>(data, size);
	}
}

MATRICE_HOST_INL decltype(auto) _Deleter(stack_alloc_tag) noexcept {
	return [](auto _dummy_) {};
}

MATRICE_HOST_INL decltype(auto) _Deleter(heap_alloc_tag) noexcept {
	/*return [](auto _ptr_) {
		if (is_aligned(_ptr_)) aligned_free(_ptr_);
		else std::free(_ptr_);
	};*/
	return [](auto _ptr_, auto _size_) {
		_Free(_ptr_, _size_);
	};
}

template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt 
_Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, stack_alloc_tag, stack_alloc_tag) {
	for (; _First != _Last; ++_Dest, (void)++_First) {
		*_Dest = *_First;
	}
	return (_Dest);
}

template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt 
_Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, stack_alloc_tag, heap_alloc_tag) {
	for (; _First != _Last; ++_Dest, (void)++_First) {
		*_Dest = *_First;
	}
	return (_Dest);
}

template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt 
_Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, heap_alloc_tag, stack_alloc_tag) {
	for (; _First != _Last; ++_Dest, (void)++_First) {
		*_Dest = *_First;
	}
	return (_Dest);
}

template<typename _InIt, typename _OutIt>
MATRICE_HOST_INL _OutIt _Memcpy(_InIt _First, _InIt _Last, _OutIt _Dest, heap_alloc_tag, heap_alloc_tag) {
	if (std::is_trivially_assignable_v<_OutIt, _InIt> || is_same_v<_InIt, _OutIt>) {
		auto const _First_ch = const_cast<const char*>(reinterpret_cast<const volatile char*>(_First));
		auto const _Last_ch = const_cast<const char*>(reinterpret_cast<const volatile char *>(_Last));
		auto const _Dest_ch = const_cast<char*>(reinterpret_cast<volatile char*>(_Dest));
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

#ifdef MATRICE_ENABLE_CUDA
MATRICE_HOST_INL decltype(auto) _Deleter(device_alloc_tag) noexcept {
	static_assert(true, "Oops, the deleter is not implemented yet.");
}

MATRICE_HOST_INL decltype(auto) _Deleter(global_alloc_tag) noexcept {
	static_assert(true, "Oops, the deleter is not implemented yet.");
}
#endif

}
}
DGE_MATRICE_END

#undef MATRICE_ALLOCTOR_SIG
#undef MATRICE_MEMCPY_ADAPTER_1
#undef MATRICE_MEMCPY_ADAPTER_2