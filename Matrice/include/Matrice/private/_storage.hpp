/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
#include <initializer_list>
#include <type_traits>
#include "../util/_macros.h"
#include "../private/_memory.h"
#include "../private/_unified_memory.h"
#include "../private/_type_traits.h"

#ifndef MATRICE_ALIGNED_CLASS
#define MATRICE_ALIGNED_CLASS class alignas(MATRICE_ALIGN_BYTES)
#endif
#ifndef MATRICE_ALIGNED_STRUCT
#define MATRICE_ALIGNED_STRUCT class alignas(MATRICE_ALIGN_BYTES)
#endif

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty> class Storage_
{
public:
#ifdef __CXX11__
	using value_t   =             _Ty;
	using pointer   =            _Ty*;
	using reference =            _Ty&;
	using int_t     =  std::ptrdiff_t;
	using idx_t     =          size_t;
#elif
	typedef _Ty               value_t;
	typedef _Ty*              pointer;
	typedef _Ty&            reference;
	typedef size_t              idx_t;
	typedef std::ptrdiff_t      int_t;
#endif
	enum Ownership { Owner = 1, Refer = 0, Proxy = -1, Empty = -2, Always = 9};
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
		enum { /*location = _Loc, */option = _Opt };
		MATRICE_GLOBAL_FINL DenseBase();
		MATRICE_GLOBAL_FINL DenseBase(int_t _rows, int_t _cols, pointer _data);
		MATRICE_GLOBAL_FINL DenseBase(int_t _rows, int_t _cols, pointer _data, initlist<value_t> _list);
		MATRICE_GLOBAL DenseBase(int_t _rows, int_t _cols);
		MATRICE_GLOBAL DenseBase(int_t _rows, int_t _cols, const value_t _val);
		MATRICE_GLOBAL DenseBase(const DenseBase& _other);
		MATRICE_GLOBAL DenseBase(DenseBase&& _other);
		template<Location _From, size_t _Option>
		MATRICE_GLOBAL_FINL DenseBase(const DenseBase<_From, _Option>& _other, pointer _data);
		template<Location _From, size_t _Option>
		MATRICE_GLOBAL_FINL DenseBase(const DenseBase<_From, _Option>& _other);
		MATRICE_GLOBAL_FINL ~DenseBase();

		///<brief> operators </brief>
		MATRICE_GLOBAL DenseBase& operator=(const DenseBase& _other);
		MATRICE_GLOBAL DenseBase& operator=(DenseBase&& _other);
		MATRICE_GLOBAL_INL decltype(auto) operator=(initlist<value_t> _list);

		///<brief> methods </brief>
		MATRICE_GLOBAL DenseBase& create(int_t _Rows, int_t _Cols);
		MATRICE_GLOBAL_FINL constexpr int_t& size() const { return (my_size);}
		MATRICE_GLOBAL_FINL constexpr int_t& rows() const { return (my_rows); }
		MATRICE_GLOBAL_FINL constexpr int_t& cols() const { return (my_cols); }
		MATRICE_GLOBAL_FINL constexpr size_t pitch() const { return (my_pitch); }
		MATRICE_GLOBAL_FINL constexpr pointer data() const { return (my_data); }
		MATRICE_GLOBAL_FINL constexpr Ownership& owner() const { return my_owner; }
		MATRICE_GLOBAL_FINL void reset(int_t rows, int_t cols, pointer data) {
			my_rows = rows, my_cols = cols, my_data = data;
		}
		MATRICE_GLOBAL_FINL constexpr bool shared() const {
#if MATRICE_SHARED_STORAGE == 1
			return (my_shared.get());
#else
			return std::false_type::value;
#endif
		}
		MATRICE_GLOBAL_FINL void free() {
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

	///<brief> generic allocator and its specialization </brief>
	template<int _M, int _N, size_t _Opt> class Allocator;

	//<brief> Managed host memory allocator </brief>
	template<int _M, int _N, size_t _Opt = allocator_traits<_M, _N>::value>
	MATRICE_ALIGNED_CLASS Allocator //: public DenseBase<OnStack, _Opt>
	{
		//typedef DenseBase<OnStack, _Opt> Base;
		using _Myt = Allocator;
	public:
		enum { location = Location::OnStack, option = _Opt };
		MATRICE_GLOBAL_INL constexpr Allocator(int ph1=0, int ph2=0) {}
			//: Base(_M, _N, _Data) {}
		MATRICE_GLOBAL_INL constexpr Allocator(int ph1, int ph2, pointer data) {
			this->_Fill_n(data);
		}
		MATRICE_GLOBAL_INL constexpr Allocator(int ph1, int ph2, value_t data) {
			this->_Fill_n(data);
		}
			//: Base(_M, _N, privt::fill_mem(data,_Data, _M*_N)) {}
		MATRICE_GLOBAL_INL constexpr Allocator(initlist<value_t> _List) {
			this->_Fill_n(_List.begin());
		}
			//: Base(_M, _N, _Data, _list) {}
		MATRICE_GLOBAL_INL constexpr Allocator(const _Myt& _other) {
			this->_Fill_n(_other._Data);
		}
			//: Base(_M, _N, privt::fill_mem(_other._Data, _Data, _other.my_size)) {}
		MATRICE_GLOBAL_INL constexpr Allocator(_Myt&& _other) {
			_Data = _other._Data;
			_other._Data = nullptr;
		}
			//: Base(std::move(_other)) {}
		//template<typename... _Args>
		//MATRICE_HOST_FINL constexpr Allocator(const _Args&... _args) 
			//: Base(_args..., _Data) {}

		MATRICE_GLOBAL_INL constexpr auto(data)() noexcept { return (_Data); }
		MATRICE_GLOBAL_INL constexpr auto(data)() const noexcept { return (_Data); }

		MATRICE_GLOBAL_INL constexpr _Myt& operator= (const _Myt& _other) noexcept {
			this->_Fill_n(_other._Data);
			/*Base::my_data = _Data;
			Base::my_cols = _other.my_cols;
			Base::my_rows = _other.my_rows;
			Base::my_size = _other.my_size;
			Base::my_owner = Owner;
			privt::fill_mem(_other._Data, _Data, Base::my_size);*/
			return (*this);
		}
		MATRICE_GLOBAL_INL constexpr _Myt& operator= (_Myt&& _other) noexcept {
			_Data = _other._Data;
			_other._Data = nullptr;
			return (*this);
		}
		MATRICE_GLOBAL_INL constexpr _Myt& operator=(const value_t _Value) noexcept {
			this->_Fill_n(_Value);
			return (*this);
		}
		MATRICE_GLOBAL_INL constexpr _Myt& operator=(const initlist<value_t> _List) noexcept {
			this->_Fill_n(_List.begin());
			return (*this);
		}
		template<typename _InIt>
		MATRICE_GLOBAL_INL constexpr _Myt& operator=(const _InIt _First) noexcept {
			this->_Fill_n(_First);
			return (*this);
		}

		MATRICE_GLOBAL_FINL constexpr auto(rows)()const noexcept { return _Myrows; }
		MATRICE_GLOBAL_FINL constexpr auto(cols)()const noexcept { return _Mycols; }
		MATRICE_GLOBAL_FINL constexpr auto(size)()const noexcept { return (_M*_N); }
		MATRICE_GLOBAL_FINL constexpr auto(pitch)()const noexcept { return one<size_t>; }
		MATRICE_GLOBAL_FINL constexpr auto(owner)()const noexcept { return (_Myown); }

	private:
		size_t _Myrows = _M, _Mycols = _N;
		Location _Myloc = Location::OnStack;
		Ownership _Myown = Ownership::Always;

		value_t _Data[_M*_N];

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
		typedef DenseBase<OnHeap, allocator_traits<0,0>::value> Base;
	public:
		enum { location = Base::location, option = Base::option };
		MATRICE_HOST_FINL Allocator(int _m, int _n) : Base(_m, _n) {}
		MATRICE_HOST_FINL Allocator(int _m, int _n, const value_t _val) : Base(_m, _n, _val) {}
		MATRICE_HOST_FINL Allocator(int _m, int _n, pointer data) : Base(_m, _n, data) {}
		MATRICE_HOST_INL Allocator(const Allocator& _other) : Base(_other) {}
		MATRICE_HOST_INL Allocator(Allocator&& _other) : Base(move(_other)) {}
		MATRICE_HOST_INL Allocator(initlist<value_t> _list) {}
		template<typename... _Args>
		MATRICE_HOST_INL Allocator(const _Args&... _args) : Base(_args...) {}

		MATRICE_HOST_INL Allocator& operator= (const initlist<value_t> _list) { return static_cast<Allocator&>(Base::operator= (_list)); }
		MATRICE_HOST_INL Allocator& operator= (const Allocator& _other) { return static_cast<Allocator&>(Base::operator= (_other)); }
		MATRICE_HOST_INL Allocator& operator= (Allocator&& _other) { return static_cast<Allocator&>(Base::operator=(move(_other))); }
	};

	//<brief> Unified device memory allocator </brief>
	template<size_t _Opt>
	MATRICE_ALIGNED_CLASS Allocator<-1, 0, _Opt> 
		: public DenseBase<OnGlobal, allocator_traits<-1, 0>::value>
	{
		typedef DenseBase<OnGlobal, allocator_traits<-1, 0>::value> Base;
	public:
		enum { location = Base::location, option = Base::option };
		MATRICE_HOST_INL Allocator() : Base() {}
		MATRICE_HOST_INL Allocator(int _m, int _n) : Base(_m, _n) {}
		MATRICE_HOST_INL Allocator(int _m, int _n, pointer data) : Base(_m, _n, data) {}
		MATRICE_HOST_INL Allocator(int _m, int _n, value_t _val) : Base(_m, _n, _val) {}
		MATRICE_HOST_INL Allocator(const Allocator& _other) : Base(_other) {}
		MATRICE_HOST_INL Allocator(Allocator&& _other) : Base(move(_other)) {}
		MATRICE_HOST_INL Allocator(initlist<value_t> _list) {}

		MATRICE_HOST_INL Allocator& operator= (const Allocator& _other) { return static_cast<Allocator&>(Base::operator= (_other)); }
	};

	//<brief> Device memory allocator </brief>
	template<size_t _Opt>
	MATRICE_ALIGNED_CLASS Allocator<-1, -1, _Opt> 
		: public DenseBase<OnDevice, allocator_traits<-1, -1>::value>
	{
		typedef DenseBase<OnDevice, allocator_traits<-1, -1>::value> Base;
	public:
		enum { location = Base::location, option = Base::option };
		MATRICE_DEVICE_INL Allocator() : Base() {}
		MATRICE_DEVICE_INL Allocator(int _m, int _n) : Base(_m, _n) {}
		MATRICE_DEVICE_INL Allocator(int _m, int _n, pointer data) : Base(_m, _n, data) {}
		MATRICE_DEVICE_INL Allocator(int _m, int _n, value_t _val) : Base(_m, _n, _val) {}
		MATRICE_DEVICE_INL Allocator(const Allocator& _other) : Base(_other) {}
		MATRICE_DEVICE_INL Allocator(Allocator&& _other) : Base(move(_other)) {}
		MATRICE_DEVICE_INL Allocator(initlist<value_t> _list) {}
		template<typename... _Args>
		MATRICE_HOST_INL Allocator(const _Args&... _args) : Base(_args...) {}

		MATRICE_DEVICE_INL Allocator& operator= (const Allocator& _other) { return static_cast<Allocator&>(Base::operator= (_other)); }
		MATRICE_GLOBAL_INL size_t pitch() const { return Base::my_pitch; }
	};
};
_DETAIL_END 
DGE_MATRICE_END
#include "inl\_storage_base.inl"