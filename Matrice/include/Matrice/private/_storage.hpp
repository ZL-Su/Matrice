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

#include <memory>
#include "private/_memory.h"
#include "private/_unified_memory.h"
#include "private/_type_traits.h"

#ifdef _MSC_VER
#pragma warning(push)
#endif

#ifndef MATRICE_ALIGNED_CLASS
#define MATRICE_ALIGNED_CLASS class alignas(MATRICE_ALIGN_BYTES)
#endif
#ifndef MATRICE_ALIGNED_STRUCT
#define MATRICE_ALIGNED_STRUCT class alignas(MATRICE_ALIGN_BYTES)
#endif

DGE_MATRICE_BEGIN
template<int _M, int _N = _M> struct allocator_traits {
	enum {
		value =
#ifdef MATRICE_ENABLE_CUDA
		(_M == ::global && _N == ::global) ? LINEAR :  // linear device allocator
		(_M == ::device && _N == ::device) ? PITCHED :  // pitched device allocator
#endif
#if MATRICE_SHARED_STORAGE == 1
		LINEAR + SHARED  // smart heap or global allocator
#else
		LINEAR + COPY    // deep heap or global allocator
#endif      
	};
};
_DETAIL_BEGIN
/**
 *\brief base class of dense allocator for plain objects.
 *\param <_Altrs> allocator traits
 */
template<class _Altrs>
MATRICE_ALIGNED_CLASS _Dense_allocator_base {
	using _Myt = _Dense_allocator_base;
	using _Mytraits = _Altrs;
public:
	using value_type = typename _Mytraits::value_type;
	using pointer = value_type*;
	using allocator = typename _Mytraits::type;
	using category = typename _Mytraits::category;
	static constexpr auto rows_at_compiletime = _Mytraits::rows;
	static constexpr auto cols_at_compiletime = _Mytraits::cols;

	MATRICE_GLOBAL_INL _Dense_allocator_base()
		: m_rows{ rows_at_compiletime }, m_cols{ cols_at_compiletime } {
		derived()._Alloc();
	}
	MATRICE_GLOBAL_INL _Dense_allocator_base(const pointer data)
		: m_rows{ rows_at_compiletime }, m_cols{ cols_at_compiletime } {
		static_assert(rows_at_compiletime*cols_at_compiletime > 0, 
			"The ctor in _Dense_allocator_base<_Altrs> is only valid for stack malloc.");
		(*this) = data;
	}
	MATRICE_GLOBAL_INL _Dense_allocator_base(size_t rows, size_t cols)
		: m_rows{ rows }, m_cols{ cols } {
		derived()._Alloc();
	}
	MATRICE_GLOBAL_INL _Dense_allocator_base(const _Myt& other)
		: _Dense_allocator_base{ other.rows(), other.cols() } {
		_Alloc_copy(other.derived());
	}
	MATRICE_GLOBAL_INL _Dense_allocator_base(_Myt&& other) noexcept {
		_Alloc_move(move(other));
	}
	template<typename _Al, enable_if_t<is_not_same_v<_Al, allocator>>>
	MATRICE_GLOBAL_INL _Dense_allocator_base(const _Al& other)
		:_Dense_allocator_base{ other.rows(), other.cols() } {
		_Alloc_copy(other);
	}
	template<typename _Al, enable_if_t<is_not_same_v<_Al, allocator>>>
	MATRICE_GLOBAL_INL _Dense_allocator_base(_Al&& other) noexcept {
		_Alloc_move(move(other));
	}

	MATRICE_GLOBAL_INL ~_Dense_allocator_base() noexcept { 
		destroy(); 
	}

	/**
	 *\brief retrieves the pointer to this memory block.
	 *\param [none]
	 */
	MATRICE_GLOBAL_INL constexpr decltype(auto)data()const noexcept {
		return (m_data);
	}
	MATRICE_GLOBAL_INL constexpr decltype(auto)data()noexcept {
		return (m_data);
	}

	MATRICE_GLOBAL_INL constexpr size_t(rows)()const noexcept {
		return m_rows;
	}
	MATRICE_GLOBAL_INL constexpr size_t(cols)()const noexcept {
		return m_cols;
	}
	MATRICE_GLOBAL_INL constexpr size_t(size)()const noexcept {
		return (m_rows*m_cols);
	}

	/**
	 *\brief returns storage order which is 101 (for row-major) or 102 (for col-major).
	 */
	MATRICE_GLOBAL_INL constexpr decltype(auto)fmt()const noexcept {
		return _Mytraits::layout_type::value;
	}

	/**
	 *\brief returns derived allocator.
	 */
	MATRICE_GLOBAL_INL const allocator& derived() const noexcept {
		return *static_cast<const allocator*>(this);
	}
	MATRICE_GLOBAL_INL allocator& derived() noexcept {
		return *static_cast<allocator*>(this);
	}

	/**
	 *\brief allocates a memory block.
	 *\param [rows, cols] rows and columns of the block to be allocated.
	 */
	MATRICE_GLOBAL_INL allocator& alloc(size_t rows=0, size_t cols=0) {
		if constexpr (is_not_same_v<category, stack_alloc_tag>) {
			m_rows = rows, m_cols = cols;
		}
		return this->derived()._Alloc();
	}

	/**
	 *\brief copy from another allocator.
	 *\param [othr] any compatible allocators.
	 */
	MATRICE_GLOBAL_INL allocator& operator=(const allocator& othr) noexcept;

	/**
	 *\brief move from another allocator.
	 *\param [othr] any compatible allocators.
	 */
	MATRICE_GLOBAL_INL allocator& operator=(allocator&& othr) noexcept;

	/**
	 *\brief fill this with a given value.
	 *\param [val] the value being filled.
	 */
	MATRICE_GLOBAL_INL allocator& operator=(const value_type val) noexcept;

	/**
	 *\brief copy from another memory block, note that the storage orders are the same and the size of source memory must not less than this->size().
	 *\param [data] the pointer to source memory
	 */
	MATRICE_GLOBAL_INL allocator& operator=(const pointer data) noexcept;

	/**
	 *\brief return if the allocator is empty or not.
	 */
	MATRICE_GLOBAL_INL operator bool() const noexcept {
		return m_data ? true : false;
	}

	/**
	 *\brief return the pointer of the allocator.
	 */
	MATRICE_GLOBAL_INL operator pointer() noexcept {
		return m_data;
	}

	/**
	 *\brief release the allocator.
	 */
	MATRICE_GLOBAL_INL void destroy() noexcept;

	/**
	 *\brief get deleter of the allocator, this method is reserved for smart pointers.
	 */
	MATRICE_GLOBAL_INL decltype(auto) deleter() const noexcept;

	/**
	 *\brief Release the ownership of the allocator
	 */
	MATRICE_GLOBAL_FINL void _Free_ownership() noexcept {
		m_data = nullptr;
		m_rows = m_cols = 0;
	}

private:
	template<typename _Al>
	MATRICE_GLOBAL_INL allocator& _Alloc_copy(const _Al& al) noexcept;
	template<typename _Al>
	MATRICE_GLOBAL_INL allocator& _Alloc_move(_Al&& al) noexcept;

protected:
	pointer m_data = nullptr;
	size_t m_rows, m_cols;
};

/**
 *\brief allocator for alloc linear stack memory.
 */
template<typename _Ty,
	diff_t _RowsAtCT, diff_t _ColsAtCT,
	size_t _Opt = allocator_traits_v<_RowsAtCT, _ColsAtCT>,
	typename _Layout = plain_layout::row_major>
MATRICE_ALIGNED_CLASS _Allocator MATRICE_NONHERITABLE
	: public _Dense_allocator_base<_Allocator_traits<_Allocator<_Ty, _RowsAtCT, _ColsAtCT, _Opt, _Layout>>> {
	using _Myt = _Allocator;
	using _Mybase = _Dense_allocator_base<_Allocator_traits<_Myt>>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::pointer;
	using _Mybase::operator=;

	MATRICE_GLOBAL_INL _Allocator() noexcept
		:_Mybase() {
	}
	MATRICE_GLOBAL_INL _Allocator(int, int) noexcept
		:_Allocator() {
	}
	template<typename _Argt>
	MATRICE_HOST_INL _Allocator(_Argt&& arg) noexcept
		:_Mybase(forward<_Argt>(arg)) {
	}

	MATRICE_HOST_INL _Myt& operator=(const _Myt& othr) noexcept {
		return _Mybase::operator=(othr);
	}

	MATRICE_GLOBAL_INL constexpr size_t(rows)()const noexcept {
		return _Mybase::rows_at_compiletime;
	}
	MATRICE_GLOBAL_INL constexpr size_t(cols)()const noexcept {
		return _Mybase::cols_at_compiletime;
	}
	MATRICE_GLOBAL_INL constexpr size_t(size)()const noexcept {
		return (rows()*cols());
	}

public:
	///</brief> reserved for internal using </brief>
	MATRICE_GLOBAL_INL constexpr decltype(auto)_Alloc() noexcept {
		_Mybase::data() = _Data; return (*this);
	}

private:
	value_type _Data[_RowsAtCT*_ColsAtCT];
};
template<typename _Ty, diff_t _M, diff_t _N, size_t _Opt, typename _Ly>
struct _Allocator_traits<_Allocator<_Ty, _M, _N, _Opt, _Ly>> {
	using value_type = _Ty;
	using layout_type = _Ly;
	using category = stack_alloc_tag;
	static constexpr auto rows = _M;
	static constexpr auto cols = _N;
	static constexpr auto options = _Opt;

	using type = detail::_Allocator<value_type, rows, cols, options, layout_type>;
};

/**
 *\brief linear memory allocator on host heap.
 */
template<typename _Ty, typename _Layout>
MATRICE_ALIGNED_CLASS _Allocator<_Ty, ::dynamic, ::dynamic, allocator_traits_v<::dynamic>, _Layout> MATRICE_NONHERITABLE
	:public _Dense_allocator_base<_Allocator_traits<_Allocator<_Ty, ::dynamic, ::dynamic, allocator_traits_v<::dynamic>, _Layout>>>
{
	using _Myt = _Allocator;
	using _Mybase = _Dense_allocator_base<_Allocator_traits<_Myt>>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::pointer;
	using _Mybase::_Dense_allocator_base;
	using _Mybase::operator=;

	MATRICE_HOST_INL _Allocator(const _Myt& other) noexcept
		:_Mybase(other) {
	}
	MATRICE_HOST_INL _Allocator(_Myt&& other) noexcept
		:_Mybase(move(other)) {
	}
	template<typename _Argt>
	MATRICE_HOST_INL _Allocator(_Argt&& arg) noexcept
		:_Mybase(forward<_Argt>(arg)){
	}
	//MATRICE_HOST_INL ~_Allocator();

	MATRICE_HOST_INL _Myt& operator=(const _Myt& othr) noexcept {
		return _Mybase::operator=(othr);
	}

public:
	MATRICE_HOST_INL decltype(auto) _Alloc() noexcept;
};
template<typename _Ty, size_t _Opt, typename _Ly>
struct _Allocator_traits<_Allocator<_Ty, ::dynamic, ::dynamic, _Opt, _Ly>> {
	using value_type = _Ty;
	using layout_type = _Ly;
	using category = heap_alloc_tag;
	static constexpr auto rows = ::dynamic;
	static constexpr auto cols = ::dynamic;
	static constexpr auto options = _Opt;

	using type = detail::_Allocator<value_type, rows, cols, options, layout_type>;
};

template<typename _Ty, diff_t _RowsAtCT, typename _Layout>
MATRICE_ALIGNED_CLASS _Allocator<_Ty, _RowsAtCT, ::dynamic, allocator_traits_v<::dynamic>, _Layout> MATRICE_NONHERITABLE
	:public _Dense_allocator_base<_Allocator_traits<_Allocator<_Ty, _RowsAtCT, ::dynamic, allocator_traits_v<::dynamic>, _Layout>>>
{
	using _Myt = _Allocator;
	using _Mybase = _Dense_allocator_base<_Allocator_traits<_Myt>>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::pointer;
	using _Mybase::_Dense_allocator_base;
	using _Mybase::operator=;

	MATRICE_HOST_INL _Allocator(const _Myt& other) noexcept
		:_Mybase(other) {
	}
	MATRICE_HOST_INL _Allocator(_Myt&& other) noexcept
		:_Mybase(move(other)) {
	}
	template<typename _Argt>
	MATRICE_HOST_INL _Allocator(_Argt&& arg) noexcept
		:_Mybase(forward<_Argt>(arg)) {
	}
	//MATRICE_HOST_INL ~_Allocator();

	MATRICE_HOST_INL _Myt& operator=(const _Myt& othr) noexcept {
		return _Mybase::operator=(othr);
	}

	MATRICE_HOST_INL constexpr size_t(rows)()const noexcept {
		return _Mybase::rows_at_compiletime;
	}

public:
	MATRICE_HOST_INL decltype(auto) _Alloc() noexcept;
};
template<typename _Ty, diff_t _RowsAtCT, size_t _Opt, typename _Ly>
struct _Allocator_traits<_Allocator<_Ty, _RowsAtCT, ::dynamic, _Opt, _Ly>> {
	using value_type = _Ty;
	using layout_type = _Ly;
	using category = heap_alloc_tag;
	static constexpr auto rows = _RowsAtCT;
	static constexpr auto cols = ::dynamic;
	static constexpr auto options = _Opt;

	using type = detail::_Allocator<value_type, rows, cols, options, layout_type>;
};

template<typename _Ty, diff_t _ColsAtCT, typename _Layout>
MATRICE_ALIGNED_CLASS _Allocator<_Ty, ::dynamic, _ColsAtCT, allocator_traits_v<::dynamic>, _Layout> MATRICE_NONHERITABLE
	:public _Dense_allocator_base<_Allocator_traits<_Allocator<_Ty, ::dynamic, _ColsAtCT, allocator_traits_v<::dynamic>, _Layout>>>
{
	using _Myt = _Allocator;
	using _Mybase = _Dense_allocator_base<_Allocator_traits<_Myt>>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::pointer;
	using _Mybase::_Dense_allocator_base;
	using _Mybase::operator=;

	MATRICE_HOST_INL _Allocator(const _Myt& other) noexcept
		:_Mybase(other) {
	}
	MATRICE_HOST_INL _Allocator(_Myt&& other) noexcept
		:_Mybase(move(other)) {
	}
	template<typename _Argt>
	MATRICE_HOST_INL _Allocator(_Argt&& arg) noexcept
		:_Mybase(forward<_Argt>(arg)) {
	}
	//MATRICE_HOST_INL ~_Allocator();

	MATRICE_HOST_INL _Myt& operator=(const _Myt& othr) noexcept {
		return _Mybase::operator=(othr);
	}

	MATRICE_HOST_INL constexpr size_t(cols)()const noexcept {
		return _Mybase::cols_at_compiletime;
	}

public:
	MATRICE_HOST_INL decltype(auto) _Alloc() noexcept;
};
template<typename _Ty, diff_t _ColsAtCT, size_t _Opt, typename _Ly>
struct _Allocator_traits<_Allocator<_Ty, ::dynamic, _ColsAtCT, _Opt, _Ly>> {
	using value_type = _Ty;
	using layout_type = _Ly;
	using category = heap_alloc_tag;
	static constexpr auto rows = ::dynamic;
	static constexpr auto cols = _ColsAtCT;
	static constexpr auto options = _Opt;

	using type = detail::_Allocator<value_type, rows, cols, options, layout_type>;
};

#ifdef MATRICE_ENABLE_CUDA
/**
 *\brief memory allocator on CUDA supprted device.
 */
template<typename _Ty, typename _Layout>
MATRICE_ALIGNED_CLASS _Allocator<_Ty, ::device, ::device, allocator_traits_v<::device>, _Layout> MATRICE_NONHERITABLE
	:public _Dense_allocator_base<_Allocator_traits<_Allocator<_Ty, ::device, ::device, allocator_traits_v<::device>, _Layout>>>
{
	using _Myt = _Allocator;
	using _Mybase = _Dense_allocator_base<_Allocator_traits<_Myt>>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::pointer;
	using _Mybase::_Dense_allocator_base;
	using _Mybase::operator=;

	MATRICE_GLOBAL_INL _Allocator(const _Myt& other) noexcept
		:_Mybase(other) {
	}
	MATRICE_GLOBAL_INL _Allocator(_Myt&& other) noexcept = delete;

	MATRICE_GLOBAL_INL ~_Allocator();

public:
	MATRICE_GLOBAL_INL decltype(auto) _Alloc() noexcept;

private:
	size_t m_pitch = _Mybase::m_cols;
};
template<typename _Ty, size_t _Opt, typename _Ly>
struct _Allocator_traits<_Allocator<_Ty, ::device, ::device, _Opt, _Ly>> {
	using value_type = _Ty;
	using layout_type = _Ly;
	using category = device_alloc_tag;
	static constexpr auto rows = ::device;
	static constexpr auto cols = ::device;
	static constexpr auto options = _Opt;

	using type = detail::_Allocator<value_type, rows, cols, options, layout_type>;
};

/**
 *\brief memory allocator on CUDA supprted device.
 */
template<typename _Ty, typename _Layout>
MATRICE_ALIGNED_CLASS _Allocator<_Ty, ::global, ::global, allocator_traits_v<::global>, _Layout> MATRICE_NONHERITABLE
	:public _Dense_allocator_base<_Allocator_traits<_Allocator<_Ty, ::global, ::global, allocator_traits_v<::global>, _Layout>>>
{
	using _Myt = _Allocator;
	using _Mybase = _Dense_allocator_base<_Allocator_traits<_Myt>>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::pointer;
	using _Mybase::_Dense_allocator_base;
	using _Mybase::operator=;

	MATRICE_GLOBAL_INL _Allocator(const _Myt& other) noexcept
		:_Mybase(other) {
	}
	MATRICE_GLOBAL_INL _Allocator(_Myt&& other) noexcept = delete;
	MATRICE_GLOBAL_INL ~_Allocator();

public:
	MATRICE_GLOBAL_INL decltype(auto) _Alloc() noexcept;
};

template<typename _Ty, size_t _Opt, typename _Ly>
struct _Allocator_traits<_Allocator<_Ty, ::global, ::global, _Opt, _Ly>> {
	using value_type = _Ty;
	using layout_type = _Ly;
	using category = global_alloc_tag;
	static constexpr auto rows = ::global;
	static constexpr auto cols = ::global;
	static constexpr auto options = _Opt;

	using type = detail::_Allocator<value_type, rows, cols, options, layout_type>;
};
#endif

template<typename _Ty> class Storage_
{
public:
#ifdef __CXX11__
	using value_t = _Ty;
	using pointer = value_t*;
	using reference = value_t&;
	using int_t = std::ptrdiff_t;
	using idx_t = size_t;
#elif
	typedef _Ty               value_t;
	typedef _Ty*              pointer;
	typedef _Ty&            reference;
	typedef size_t              idx_t;
	typedef std::ptrdiff_t      int_t;
#endif
	enum Ownership { Owner = 1, Refer = 0, Proxy = -1, Empty = -2, Always = 9};

	///<brief> generic allocator and its specialization </brief>
	template<int _M, int _N, size_t _Opt> class Allocator;

	///<brief> data memory </brief>
	template<
		Location _Loc = UnSpecified, 
		size_t _Opt = LINEAR+
#if MATRICE_SHARED_STORAGE == 1
		SHARED
#else
		COPY
#endif 
			 >
	MATRICE_ALIGNED_CLASS DenseBase {
	public:
		static constexpr auto location = _Loc;
		enum { option = _Opt };
		MATRICE_GLOBAL_FINL DenseBase() noexcept;
		MATRICE_GLOBAL_FINL DenseBase(int_t _rows, int_t _cols, pointer _data);
		MATRICE_GLOBAL_FINL DenseBase(int_t _rows, int_t _cols, pointer _data, initlist<value_t> _list);
		MATRICE_GLOBAL DenseBase(int_t _rows, int_t _cols);
		MATRICE_GLOBAL DenseBase(int_t _rows, int_t _cols, const value_t _val);
		MATRICE_GLOBAL_FINL DenseBase(initlist<value_t> _list);
		MATRICE_GLOBAL DenseBase(const DenseBase& _other);
		MATRICE_GLOBAL DenseBase(DenseBase&& _other) noexcept;

		template<Location _From, size_t _Option>
		MATRICE_GLOBAL_FINL DenseBase(const DenseBase<_From, _Option>& _other, pointer _data);
		template<Location _From, size_t _Option>
		MATRICE_GLOBAL_FINL DenseBase(const DenseBase<_From, _Option>& _other);
		MATRICE_GLOBAL_FINL ~DenseBase() noexcept;

		///<brief> operators </brief>
		MATRICE_GLOBAL DenseBase& operator=(const DenseBase& _other);
		MATRICE_GLOBAL DenseBase& operator=(DenseBase&& _other) noexcept;
		MATRICE_GLOBAL_INL decltype(auto) operator=(value_t _value) noexcept;
		MATRICE_GLOBAL_INL decltype(auto) operator=(initlist<value_t> _list) noexcept;

		template<int _M, int _N>
		MATRICE_GLOBAL_INL DenseBase(const Allocator<_M, _N, allocator_traits<_M, _N>::value>& _al) noexcept;
		template<int _M, int _N>
		MATRICE_GLOBAL_INL decltype(auto) operator=(const Allocator<_M, _N, allocator_traits<_M, _N>::value>& _al) noexcept;

		///<brief> methods </brief>
		MATRICE_GLOBAL DenseBase& create(int_t _Rows, int_t _Cols);
		MATRICE_GLOBAL_FINL constexpr int_t& size() const noexcept { return (my_size);}
		MATRICE_GLOBAL_FINL constexpr int_t& rows() const noexcept { return (my_rows); }
		MATRICE_GLOBAL_FINL constexpr int_t& cols() const noexcept { return (my_cols); }
		MATRICE_GLOBAL_FINL constexpr size_t pitch() const noexcept { return (my_pitch); }
		MATRICE_GLOBAL_FINL constexpr pointer data() const noexcept { return (my_data); }
		MATRICE_GLOBAL_FINL constexpr Ownership& owner() const noexcept { return my_owner; }
		MATRICE_GLOBAL_FINL void reset(int_t rows, int_t cols, pointer data) noexcept {
			my_rows = rows, my_cols = cols, my_data = data;
		}
		MATRICE_GLOBAL_FINL constexpr bool shared() const noexcept {
#if MATRICE_SHARED_STORAGE == 1
			return (my_shared.get());
#else
			return std::false_type::value;
#endif
		}
		MATRICE_GLOBAL_FINL void free() noexcept {
			my_cols = 0, my_rows = 0, my_size = 0;
			my_owner = Empty, my_pitch = 0;
			my_location = UnSpecified;
#if MATRICE_SHARED_STORAGE == 1
			if (my_shared) my_shared = nullptr;
#endif
			my_data = nullptr;
		}

	protected:
		mutable int_t my_rows, my_cols;
		mutable int_t my_size;
		mutable pointer my_data;
		mutable Ownership my_owner = Empty;
		size_t my_pitch = 1; //used for CUDA pitched malloc only
	private:
#if MATRICE_SHARED_STORAGE == 1
		using SharedPtr = std::shared_ptr<value_t>;
		SharedPtr my_shared;
#endif
		Location my_location = _Loc;
	};

	//<brief> Managed host memory allocator </brief>
	template<int _M, int _N=_M, size_t _Opt = allocator_traits_v<_M, _N>>
	MATRICE_ALIGNED_CLASS Allocator {
		using _Myt = Allocator;
	public:
		enum { location = Location::OnStack, option = _Opt };

		MATRICE_GLOBAL_INL constexpr Allocator(int ph1=0, int ph2=0) noexcept {
		}
		MATRICE_GLOBAL_INL constexpr Allocator(int ph1, int ph2, pointer data) noexcept {
			this->_Fill_n(data);
		}
		MATRICE_GLOBAL_INL constexpr Allocator(int ph1, int ph2, value_t data) noexcept {
			this->_Fill_n(data);
		}
		MATRICE_GLOBAL_INL constexpr Allocator(initlist<value_t> _List) noexcept {
			this->_Fill_n(_List.begin());
		}
		MATRICE_GLOBAL_INL constexpr Allocator(const _Myt& _other) noexcept {
			this->_Fill_n(_other._Data);
		}
		MATRICE_GLOBAL_INL constexpr Allocator(_Myt&& _other) noexcept {
			this->_Fill_n(_other._Data);
		}
		template<typename _Alty>
		MATRICE_HOST_FINL constexpr Allocator(const _Alty& _alloc) noexcept {
			this->_Fill_n(_alloc.data());
		}

		MATRICE_GLOBAL_INL constexpr auto(data)() noexcept { return (_Data); }
		MATRICE_GLOBAL_INL constexpr auto(data)() const noexcept { return (_Data); }

		MATRICE_GLOBAL_INL constexpr _Myt& operator= (const _Myt& _other) noexcept {
			this->_Fill_n(_other._Data);
			return (*this);
		}
		MATRICE_GLOBAL_INL constexpr _Myt& operator= (_Myt&& _other) noexcept {
			this->_Fill_n(_other._Data);
			return (*this);
		}
		MATRICE_GLOBAL_INL constexpr _Myt& operator=(const value_t _Value) noexcept {
			this->_Fill_n(_Value);
			return (*this);
		}
		template<typename _InIt>
		MATRICE_GLOBAL_INL constexpr _Myt& operator=(const _InIt _First) noexcept {
			this->_Fill_n(_First);
			return (*this);
		}
		MATRICE_GLOBAL_INL constexpr _Myt& operator=(const initlist<value_t> _List) noexcept {
			if (_List.size() == 1) this->operator=(*_List.begin());
			else this->operator=(_List.begin());
			return (*this);
		}

		MATRICE_GLOBAL_FINL constexpr auto(rows)()const noexcept { return _Myrows; }
		MATRICE_GLOBAL_FINL constexpr auto(cols)()const noexcept { return _Mycols; }
		MATRICE_GLOBAL_FINL constexpr auto(size)()const noexcept { return (_M*_N); }
		MATRICE_GLOBAL_FINL constexpr auto(owner)()const noexcept { return (_Myown); }
		MATRICE_GLOBAL_FINL constexpr auto(pitch)()const noexcept { return size_t(1); }

		// \fill _Data with zeros, because its managed.
		MATRICE_GLOBAL_FINL constexpr void free() noexcept {
			this->_Fill_n(value_t(0));
		}

		// \reserved method for using smart pointer.
		MATRICE_GLOBAL_FINL constexpr auto deleter() noexcept {
			return [](auto _) {};
		}

	private:
		size_t _Myrows = _M, _Mycols = _N;
		Location _Myloc = Location::OnStack;
		Ownership _Myown = Ownership::Always;

		value_t _Data[_M*_N];

		/**
		 *\brief fill managed host memory.
		 *\param [_First] can be a scalar or an iterator/pointer type.
		 */
		template<typename _InIt>
		MATRICE_GLOBAL_INL constexpr void _Fill_n(_InIt _First) noexcept {
			for (auto _Idx = 0; _Idx < size(); ++_Idx) {
				if constexpr (is_iterator_v<_InIt>)
					_Data[_Idx] = _First[_Idx];
				else
					_Data[_Idx] = _First;
			}
		}
	};

	//<brief> Dynamic host memory allocator </brief>
	template<size_t _Opt>
	MATRICE_ALIGNED_CLASS Allocator<::dynamic, ::dynamic, _Opt>
		: public DenseBase<OnHeap, _Opt>
	{
		using _Mybase = DenseBase<OnHeap, allocator_traits<::dynamic>::value>;
	public:
		enum { location = _Mybase::location, option = _Mybase::option };
		MATRICE_HOST_FINL Allocator(int _m, int _n) 
			:_Mybase(_m, _n) {}
		MATRICE_HOST_FINL Allocator(int _m, int _n, const value_t _val) 
			:_Mybase(_m, _n, _val) {}
		MATRICE_HOST_FINL Allocator(int _m, int _n, pointer data) 
			:_Mybase(_m, _n, data) {}
		MATRICE_HOST_INL Allocator(const Allocator& _other) 
			:_Mybase(_other) {}
		MATRICE_HOST_INL Allocator(Allocator&& _other) 
			:_Mybase(move(_other)) {}
		MATRICE_HOST_INL Allocator(initlist<value_t> _list)
			:_Mybase(_list) {}
		template<typename... _Args>
		MATRICE_HOST_INL Allocator(_Args&&... _args) 
			:_Mybase(forward<_Args>(_args)...) {}

		template<typename... _Args>
		MATRICE_HOST_INL Allocator& operator=(_Args&&... args) noexcept {
			return static_cast<Allocator&>(_Mybase::operator=(forward<_Args>(args)...));
		}
		MATRICE_HOST_INL Allocator& operator= (const Allocator& _other) { 
			return static_cast<Allocator&>(_Mybase::operator= (_other)); 
		}
		MATRICE_HOST_INL Allocator& operator= (Allocator&& _other) { 
			return static_cast<Allocator&>(_Mybase::operator=(move(_other))); 
		}

		// \reserved method for using smart pointer.
		MATRICE_GLOBAL_FINL constexpr auto deleter() noexcept {
				return [](auto _) { 
					if (internal::_Memory::is_aligned(_))
						privt::aligned_free(_);
					else
						std::free(_);
				};
		}
	};

#ifdef MATRICE_ENABLE_CUDA
	//<brief> Unified device memory allocator </brief>
	template<size_t _Opt>
	MATRICE_ALIGNED_CLASS Allocator<::global, ::global, _Opt> 
		: public DenseBase<OnGlobal, allocator_traits<::global>::value>
	{
		using _Mybase = DenseBase<OnGlobal, allocator_traits<::global>::value>;
	public:
		enum { location = _Mybase::location, option = _Mybase::option };
		MATRICE_HOST_INL Allocator() 
			: _Mybase() {}
		MATRICE_HOST_INL Allocator(int _m, int _n) 
			: _Mybase(_m, _n) {}
		MATRICE_HOST_INL Allocator(int _m, int _n, pointer data) 
			: _Mybase(_m, _n, data) {}
		MATRICE_HOST_INL Allocator(int _m, int _n, value_t _val) 
			: _Mybase(_m, _n, _val) {}
		MATRICE_HOST_INL Allocator(const Allocator& _other) 
			: _Mybase(_other) {}
		MATRICE_HOST_INL Allocator(Allocator&& _other) 
			: _Mybase(move(_other)) {}
		MATRICE_HOST_INL Allocator(initlist<value_t> _list) 
			: _Mybase(_list) {}
		template<typename... _Args>
		MATRICE_HOST_INL Allocator(_Args&&... _args)
			: _Mybase(forward<_Args>(_args)...) {}

		MATRICE_HOST_INL Allocator& operator= (const Allocator& _other) { 
			return static_cast<Allocator&>(_Mybase::operator= (_other)); 
		}

		// \reserved method for using smart pointer.
		MATRICE_GLOBAL_INL constexpr auto deleter() noexcept {
			return [](auto _) { privt::device_free(_); };
		}
	};

	//<brief> Device memory allocator </brief>
	template<size_t _Opt>
	MATRICE_ALIGNED_CLASS Allocator<::device, ::device, _Opt> 
		: public DenseBase<OnDevice, allocator_traits<::device>::value>
	{
		using _Mybase = DenseBase<OnDevice, allocator_traits<::device>::value>;
	public:
		enum { location = _Mybase::location, option = _Mybase::option };
		MATRICE_DEVICE_INL Allocator() 
			: _Mybase() {}
		MATRICE_DEVICE_INL Allocator(int _m, int _n) 
			: _Mybase(_m, _n) {}
		MATRICE_DEVICE_INL Allocator(int _m, int _n, pointer data) 
			: _Mybase(_m, _n, data) {}
		MATRICE_DEVICE_INL Allocator(int _m, int _n, value_t _val) 
			: _Mybase(_m, _n, _val) {}
		MATRICE_DEVICE_INL Allocator(const Allocator& _other) 
			: _Mybase(_other) {}
		MATRICE_DEVICE_INL Allocator(Allocator&& _other) 
			: _Mybase(move(_other)) {}
		MATRICE_DEVICE_INL Allocator(initlist<value_t> _list)
			: _Mybase(_list) {}
		template<typename... _Args>
		MATRICE_HOST_INL Allocator(_Args&&... _args)
			: _Mybase(forward<_Args>(_args)...) {}

		MATRICE_DEVICE_INL Allocator& operator=(const Allocator& _other) { 
			return static_cast<Allocator&>(_Mybase::operator= (_other)); 
		}
		MATRICE_GLOBAL_INL size_t pitch() const { 
			return _Mybase::my_pitch; 
		}

		// \reserved method for using smart pointer.
		MATRICE_GLOBAL_INL constexpr auto deleter() noexcept {
			return [](auto _) { privt::device_free(_); };
		}
	};
#endif 
};
_DETAIL_END

/**
 *\brief dense allocator interface
 *\param <see comments below>
 *\note if both _RowsAtCT and _ColsAtCT are greater than zero, the allocator will create memory block on the stack; if they are zero, the allocator will be dynamically create memory block on the heap; else if they are equal to -1 or -2, the memory will be allocated on a device (GPU), or with the CUDA managed malloc technique.
 */
template<typename _Ty, /*data type*/
	diff_t _RowsAtCT = 0, /*rows at compile-time*/
	diff_t _ColsAtCT = _RowsAtCT, /*cols at compile-time*/
	class _Ly=plain_layout::row_major /*storage order*/
>
using dense_allocator = detail::_Allocator<_Ty, _RowsAtCT, _ColsAtCT, allocator_traits_v<_RowsAtCT, _ColsAtCT>, _Ly>;

namespace internal {
template<typename _Ty>
using _Dynamic_buffer = detail::_Allocator<_Ty, ::dynamic, ::dynamic, allocator_traits_v<::dynamic, ::dynamic>, plain_layout::row_major>;
}

DGE_MATRICE_END

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include "storage/_allocator.inl"
#include "inl\_storage_base.inl"