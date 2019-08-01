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
#include "../util/_macros.h"
#include "../private/_memory.h"
#include "../private/_unified_memory.h"
#include "../private/_type_traits.h"

#ifdef _MSC_VER
#pragma warning(push)
#endif

#ifndef MATRICE_ALIGNED_CLASS
#define MATRICE_ALIGNED_CLASS class alignas(MATRICE_ALIGN_BYTES)
#endif
#ifndef MATRICE_ALIGNED_STRUCT
#define MATRICE_ALIGNED_STRUCT class alignas(MATRICE_ALIGN_BYTES)
#endif
namespace {
	static constexpr int dynamic = 0;
	static constexpr int device = -1;
	static constexpr int global = -2;
}

DGE_MATRICE_BEGIN
template<int _M, int _N> struct allocator_traits {
	enum {
		value = (_M == 0 && _N == -1) ? LINEAR :  // linear device allocator
		(_M == -1 && _N == -1) ? PITCHED :  // pitched device allocator
#if MATRICE_SHARED_STORAGE == 1
		LINEAR + SHARED  // smart heap or global allocator
#else
		LINEAR + COPY    // deep heap or global allocator
#endif      
	};
};
struct allocator_tag {
	struct stack_allocator {};
	struct heap_allocator {};
	struct device_allocator {};
	struct global_allocator {};
};
struct plain_layout {
	struct row_major { static constexpr auto value = 101; };
	struct col_major { static constexpr auto value = 102; };
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
	using pointer = std::add_pointer_t<value_type>;
	using allocator = typename _Mytraits::allocator;
	static constexpr auto rows_at_compiletime = _Mytraits::rows;
	static constexpr auto cols_at_compiletime = _Mytraits::cols;

	MATRICE_GLOBAL_INL _Dense_allocator_base()
		:m_rows(rows_at_compiletime), m_cols(cols_at_compiletime){
		derived()._Alloc();
	}
	MATRICE_GLOBAL_INL _Dense_allocator_base(size_t rows, size_t cols)
		:m_rows(rows), m_cols(cols) {
		derived()._Alloc();
	}
	MATRICE_GLOBAL_INL _Dense_allocator_base(const _Myt& othr)
		: m_rows(othr.m_rows), m_cols(othr.m_cols) {
		derived()._Alloc()._Copy(othr);
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
	 \param [none]
	 */
	MATRICE_GLOBAL_INL constexpr decltype(auto)fmt()const noexcept {
		return _Mytraits::layout_type::value;
	}

	/**
	 *\brief returns derived allocator.
	 *\param [none]
	 */
	MATRICE_GLOBAL_INL const allocator& derived() const noexcept {
		return *static_cast<const allocator*>(this);
	}
	MATRICE_GLOBAL_INL allocator& derived() noexcept {
		return *static_cast<allocator*>(this);
	}

	MATRICE_GLOBAL_INL _Myt& operator=(const value_type val) noexcept {
		if (m_data != nullptr) {
			for (auto idx = 0; idx < this->size(); ++idx) {
				m_data[idx] = val;
			}
		}
		return (*this);
	}
	/**
	 *\brief copy from another memory block, note that the storage orders are the same and the size of source memory must not less than this->size().
	 *\param [data] the pointer to source memory
	 */
	MATRICE_GLOBAL_INL _Myt& operator=(const pointer data) noexcept {
		if (m_data != nullptr && data) {
			for (auto idx = 0; idx < this->size(); ++idx) {
				m_data[idx] = data[idx];
			}
		}
		return (*this);
	}
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

	MATRICE_GLOBAL_INL constexpr size_t(rows)() noexcept {
		return _Mybase::rows_at_compiletime;
	}
	MATRICE_GLOBAL_INL constexpr size_t(cols)() noexcept {
		return _Mybase::cols_at_compiletime;
	}
	MATRICE_GLOBAL_INL constexpr size_t(size)() noexcept {
		return (rows()*cols());
	}

private:
	value_type _Data[_RowsAtCT*_ColsAtCT];

public:
	///</brief> reserved for internal using </brief>
	MATRICE_GLOBAL_INL constexpr decltype(auto)_Alloc() noexcept {
		_Mybase::data() = _Data; return (*this);
	}
};

template<typename _Ty, diff_t _M, diff_t _N, 
	size_t _Opt, typename _Ly>
struct _Allocator_traits<_Allocator<_Ty, _M, _N, _Opt, _Ly>> {
	using value_type = _Ty;
	using layout_type = _Ly;
	using category = allocator_tag::stack_allocator;
	static constexpr auto rows = _M;
	static constexpr auto cols = _N;
	static constexpr auto options = _Opt;

	using allocator = detail::_Allocator<value_type, rows, cols, options, layout_type>;
};

/**
 *\brief linear memory allocator on host heap.
 */
template<typename _Ty, typename _Layout>
MATRICE_ALIGNED_CLASS _Allocator<_Ty, 0, 0, allocator_traits_v<0, 0>, _Layout> MATRICE_NONHERITABLE
	: public _Dense_allocator_base<_Allocator_traits<_Allocator<_Ty, 0, 0, allocator_traits_v<0, 0>, _Layout>>>{
	using _Myt = _Allocator;
	using _Mybase = _Dense_allocator_base<_Allocator_traits<_Myt>>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::pointer;
	using _Mybase::operator=;

	MATRICE_HOST_INL _Allocator() noexcept
		:_Mybase() {
	}
	MATRICE_HOST_INL _Allocator(size_t rows, size_t cols = 1) noexcept
		:_Mybase(rows, cols) {
	}
	MATRICE_HOST_INL ~_Allocator();

public:
	MATRICE_HOST_INL decltype(auto) _Alloc() noexcept;
	MATRICE_HOST_INL decltype(auto) _Copy(const _Mybase& othr) noexcept;
	MATRICE_HOST_INL decltype(auto) _Move(_Mybase&& othr) noexcept;
};

template<typename _Ty, size_t _Opt, typename _Ly>
struct _Allocator_traits<_Allocator<_Ty, 0, 0, _Opt, _Ly>> {
	using value_type = _Ty;
	using layout_type = _Ly;
	using category = allocator_tag::heap_allocator;
	static constexpr auto rows = 0;
	static constexpr auto cols = 0;
	static constexpr auto options = _Opt;

	using allocator = detail::_Allocator<value_type, rows, cols, options, layout_type>;
};

template<typename _Ty> class Storage_
{
public:
#ifdef __CXX11__
	using value_t = _Ty;
	using pointer = std::add_pointer_t<value_t>;
	using reference = std::add_lvalue_reference_t<value_t>;
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
	template<int _M, int _N, size_t _Opt = allocator_traits<_M, _N>::value>
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
	MATRICE_ALIGNED_CLASS Allocator<0, 0, _Opt> 
		: public DenseBase<OnHeap, _Opt>
	{
		using _Mybase = DenseBase<OnHeap, allocator_traits<0,0>::value>;
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

	//<brief> Unified device memory allocator </brief>
	template<size_t _Opt>
	MATRICE_ALIGNED_CLASS Allocator<-1, 0, _Opt> 
		: public DenseBase<OnGlobal, allocator_traits<-1, 0>::value>
	{
		using _Mybase = DenseBase<OnGlobal, allocator_traits<-1, 0>::value>;
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
	MATRICE_ALIGNED_CLASS Allocator<-1, -1, _Opt> 
		: public DenseBase<OnDevice, allocator_traits<-1, -1>::value>
	{
		using _Mybase = DenseBase<OnDevice, allocator_traits<-1, -1>::value>;
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
DGE_MATRICE_END

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include "inl\_storage_base.inl"
#include "storage\_allocator.inl"