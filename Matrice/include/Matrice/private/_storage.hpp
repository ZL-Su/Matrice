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
more details.

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
#include "../private/_expr_type_traits.h"

#ifndef MATRICE_ALIGNED_CLASS
#define MATRICE_ALIGNED_CLASS class alignas(MATRICE_ALIGN_BYTES)
#endif
#ifndef MATRICE_ALIGNED_STRUCT
#define MATRICE_ALIGNED_STRUCT class alignas(MATRICE_ALIGN_BYTES)
#endif

namespace dgelom { namespace details {
template<typename _Ty> class Storage_
{
public:
#ifdef __CXX11__
	using value_t   =             _Ty;
	using pointer   =            _Ty*;
	using reference =            _Ty&;
	using int_t     =  std::ptrdiff_t;
	using idx_t     =     std::size_t;
#elif
	typedef _Ty               value_t;
	typedef _Ty*              pointer;
	typedef _Ty&            reference;
	typedef std::size_t         idx_t;
	typedef std::ptrdiff_t      int_t;
#endif
	enum Ownership { Owner = 1, Refer = 0, Proxy = -1, Dummy = -2 };
	///<brief> data memory </brief>
	template<Location _Loc = UnSpecified, size_t _Opt = LINEAR+
#ifdef __CXX11_SHARED__
		SHARED
#else
		COPY
#endif
	> 
	MATRICE_ALIGNED_CLASS DenseBase {
	public:
		enum { location = _Loc, option = _Opt };
		MATRICE_GLOBAL_FINL DenseBase();
		MATRICE_GLOBAL_FINL DenseBase(int_t _rows, int_t _cols, pointer _data);
		MATRICE_GLOBAL_FINL DenseBase(int_t _rows, int_t _cols, pointer _data, std::initializer_list<value_t> _list);
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
		MATRICE_GLOBAL_FINL auto& operator=(std::initializer_list<value_t> _list);

		///<brief> methods </brief>
		MATRICE_GLOBAL_FINL int_t& size() const { return (my_size);}
		MATRICE_GLOBAL_FINL int_t& rows() const { return (my_rows); }
		MATRICE_GLOBAL_FINL int_t& cols() const { return (my_cols); }
		MATRICE_GLOBAL_FINL size_t pitch() const { return (my_pitch); }
		MATRICE_GLOBAL_FINL pointer data() const { return (my_data); }
		MATRICE_GLOBAL_FINL Ownership& owner() const { return my_owner; }
		MATRICE_GLOBAL_FINL bool shared() const 
		{
#ifdef __CXX11_SHARED__
			return (my_shared.get());
#else
			return 0;
#endif
		}
		MATRICE_GLOBAL_FINL void free()
		{
			my_cols = 0, my_rows = 0, my_size = 0;
			my_owner = Dummy, my_pitch = 0;
			my_location = UnSpecified;
			my_shared = nullptr;
			my_data = nullptr;
		}

	protected:
		mutable int_t my_rows, my_cols;
		mutable int_t my_size;
		mutable pointer my_data;
		mutable Ownership my_owner = Owner;
		std::size_t my_pitch = 1; //used for CUDA pitched malloc only
	private:
#ifdef __CXX11_SHARED__
		using SharedPtr = std::shared_ptr<value_t>;
		SharedPtr my_shared;
#endif
		Location my_location = _Loc;
	};

	///<brief> generic allocator and its specialization </brief>
	template<int _M, int _N, size_t _Opt> class Allocator;

	//<brief> Managed host memory allocator </brief>
	template<int _M, int _N, size_t _Opt = LINEAR+COPY>
	MATRICE_ALIGNED_CLASS Allocator : public DenseBase<OnStack, _Opt>
	{
		typedef DenseBase<OnStack, _Opt> Base;
	public:
		enum { location = Base::location, option = Base::option };
		MATRICE_HOST_FINL Allocator(int ph1=0, int ph2=0) : Base(_M, _N, _Data) {}
		MATRICE_HOST_FINL Allocator(int ph1, int ph2, pointer data) : Base(_M, _N, data) {}
		MATRICE_HOST_FINL Allocator(std::initializer_list<value_t> _list) : Base(_M, _N, _Data, _list) {}
		MATRICE_HOST_FINL Allocator(const Allocator& _other) : Base(_M, _N, privt::fill_mem(_other._Data, _Data, _other.my_size)) {}
		MATRICE_HOST_FINL Allocator(Allocator&& _other) : Base(std::move(_other)) {}
		template<typename... _Args>
		MATRICE_HOST_FINL Allocator(const _Args&... _args) : Base(_args..., _Data) {}
		MATRICE_HOST_FINL Allocator& operator= (const Allocator& _other)
		{
			Base::my_data = _Data;
			Base::my_cols = _other.my_cols;
			Base::my_rows = _other.my_rows;
			Base::my_size = _other.my_size;
			Base::my_owner = _other.my_owner;
			privt::fill_mem(_other._Data, _Data, Base::my_size);
			return (*this);
		}
	private:
		value_t _Data[_M*_N];
	};
	//<brief> Dynamic host memory allocator </brief>
	template<size_t _Opt> 
	MATRICE_ALIGNED_CLASS Allocator<0, 0, _Opt> : public DenseBase<OnHeap, _Opt>
	{
		typedef DenseBase<OnHeap, _Opt> Base;
	public:
		enum { location = Base::location, option = Base::option };
		MATRICE_HOST_INL Allocator() : Base() {}
		MATRICE_HOST_INL Allocator(int _m, int _n) : Base(_m, _n) {}
		MATRICE_HOST_INL Allocator(int _m, int _n, const value_t _val) : Base(_m, _n, _val) {}
		MATRICE_HOST_INL Allocator(int _m, int _n, pointer data) : Base(_m, _n, data) {}
		MATRICE_HOST_INL Allocator(const Allocator& _other) : Base(_other) {}
		MATRICE_HOST_INL Allocator(Allocator&& _other) : Base(std::move(_other)) {}
		MATRICE_HOST_INL Allocator(std::initializer_list<value_t> _list) {}
		template<typename... _Args>
		MATRICE_HOST_INL Allocator(const _Args&... _args) : Base(_args...) {}

		MATRICE_HOST_INL Allocator& operator= (const std::initializer_list<value_t> _list) { return static_cast<Allocator&>(Base::operator= (_list)); }
		MATRICE_HOST_INL Allocator& operator= (const Allocator& _other) { return static_cast<Allocator&>(Base::operator= (_other)); }
		MATRICE_HOST_INL Allocator& operator= (Allocator&& _other) { return static_cast<Allocator&>(Base::operator=(std::move(_other))); }
	};
	//<brief> Unified device memory allocator </brief>
	template<size_t _Opt>
	MATRICE_ALIGNED_CLASS Allocator<-1, 0, _Opt> : public DenseBase<OnGlobal, _Opt>
	{
		typedef DenseBase<OnGlobal, _Opt> Base;
	public:
		enum { location = Base::location, option = Base::option };
		MATRICE_HOST_INL Allocator() : Base() {}
		MATRICE_HOST_INL Allocator(int _m, int _n) : Base(_m, _n) {}
		MATRICE_HOST_INL Allocator(int _m, int _n, pointer data) : Base(_m, _n, data) {}
		MATRICE_HOST_INL Allocator(int _m, int _n, value_t _val) : Base(_m, _n, _val) {}
		MATRICE_HOST_INL Allocator(const Allocator& _other) : Base(_other) {}
		MATRICE_HOST_INL Allocator(Allocator&& _other) : Base(std::move(_other)) {}
		MATRICE_HOST_INL Allocator(std::initializer_list<value_t> _list) {}

		MATRICE_HOST_INL Allocator& operator= (const Allocator& _other) { return static_cast<Allocator&>(Base::operator= (_other)); }
	};
	//<brief> Device memory allocator </brief>
	template<size_t _Opt>
	MATRICE_ALIGNED_CLASS Allocator<-1, -1, _Opt> : public DenseBase<OnDevice, _Opt>
	{
		typedef DenseBase<OnDevice, _Opt> Base;
	public:
		enum { location = Base::location, option = Base::option };
		MATRICE_DEVICE_INL Allocator() : Base() {}
		MATRICE_DEVICE_INL Allocator(int _m, int _n) : Base(_m, _n) {}
		MATRICE_DEVICE_INL Allocator(int _m, int _n, pointer data) : Base(_m, _n, data) {}
		MATRICE_DEVICE_INL Allocator(int _m, int _n, value_t _val) : Base(_m, _n, _val) {}
		MATRICE_DEVICE_INL Allocator(const Allocator& _other) : Base(_other) {}
		MATRICE_DEVICE_INL Allocator(Allocator&& _other) : Base(std::move(_other)) {}
		MATRICE_DEVICE_INL Allocator(std::initializer_list<value_t> _list) {}
		template<typename... _Args>
		MATRICE_HOST_INL Allocator(const _Args&... _args) : Base(_args...) {}

		MATRICE_DEVICE_INL Allocator& operator= (const Allocator& _other) { return static_cast<Allocator&>(Base::operator= (_other)); }
		MATRICE_GLOBAL_INL std::size_t pitch() const { return Base::my_pitch; }
	};

#pragma region <!-- decrepated -->
	///<brief> specialization allocator classes </brief>
	template<class DerivedAllocator, Location _Loc> class Base_
	{
#ifdef __CXX11__
		using SharedPtr = std::shared_ptr<value_t>;
#endif
	public:
		Base_() : m_rows(0), m_cols(0), m_location(UnSpecified){}
		Base_(int rows, int cols) : m_rows(rows), m_cols(cols), m_data(privt::aligned_malloc<value_t>(rows*cols)) {}
#ifdef __CXX11__
		Base_(int rows, int cols, bool is_shared) : m_rows(rows), m_cols(cols), m_shared(is_shared)
		{
			_SharedData = SharedPtr(privt::aligned_malloc<value_t>(rows*cols), [=](pointer ptr) { privt::aligned_free(ptr); });
			m_data = _SharedData.get();
		}
#endif
		~Base_()
		{
			if ((!m_shared) && m_data) privt::aligned_free(m_data);
		}

#pragma region operator overload
		MATRICE_GLOBAL_INL
		reference operator[](idx_t idx) const { return m_data[idx]; }
#pragma endregion
#pragma region public methods
		MATRICE_GLOBAL_INL std::size_t cols() const { return m_cols; }
		MATRICE_GLOBAL_INL std::size_t rows() const { return m_rows; }
		MATRICE_GLOBAL_INL std::size_t size() const { return m_rows * m_cols; }
		MATRICE_GLOBAL_INL pointer data() const { return (m_data); }
		MATRICE_GLOBAL_INL pointer ptr(idx_t ridx) const { return (m_data + m_cols*ridx); }
		MATRICE_GLOBAL_INL loctn_t& location() const { return m_location; }
#pragma endregion

	protected:
#ifdef __CXX11__
		
		pointer m_data = nullptr;
#elif       
		pointer m_data;
#endif
		mutable Location m_location = _Loc;
		std::size_t m_rows, m_cols;
		bool m_shared = false;
	private:
#ifdef __CXX11__
		SharedPtr _SharedData = nullptr;
#endif
	};
	template<int M, int N> MATRICE_ALIGNED_CLASS ManagedAllocator : public Base_<ManagedAllocator<M, N>, OnStack>
	{
		using _Base = Base_<ManagedAllocator<M, N>, OnStack>;
		using _Base::m_rows;
		using _Base::m_cols;
		typedef pointer MyPtr;
	public:
		ManagedAllocator() : _Base(M, N) { m_data = &_Data[0]; }
		ManagedAllocator(int _m, int _n) : _Base(M, N) { m_data = &_Data[0]; }

	private:
		value_t _Data[M*N];
		using _Base::m_data;
	};
	template<typename = std::enable_if<std::is_arithmetic<_Ty>::value, _Ty>::type>
	class DynamicAllocator : public Base_<DynamicAllocator<_Ty>, OnHeap>
	{
		using _Base = Base_<DynamicAllocator<_Ty>, OnHeap>;
		using _Base::m_rows;
		using _Base::m_cols;
		typedef pointer MyPtr;
	public:
		DynamicAllocator() :_Base() {}
		template<typename IntegerType>
		DynamicAllocator(const IntegerType& _M, const IntegerType& _N) :_Base(_M, _N) {}

	private:
		using _Base::m_data;
	};
#ifdef __CXX11__
	template<typename = std::enable_if<std::is_arithmetic<_Ty>::value, _Ty>::type> 
	class SharedAllocator : public Base_<SharedAllocator<_Ty>, OnHeap>
	{
		using _Base = Base_<SharedAllocator<_Ty>, OnHeap>;
		using _Base::m_rows;
		using _Base::m_cols;
	public:
		SharedAllocator() : _Base() {}
		template<typename IntegerType>
		SharedAllocator(const IntegerType& _M, const IntegerType& _N) :_Base(_M, _N, true) {}

	private:
		using _Base::m_data;
	};
#endif
	template<typename = std::enable_if<std::is_arithmetic<_Ty>::value, _Ty>::type>
	class UnifiedAllocator : public Base_<UnifiedAllocator<_Ty>, OnGlobal>
	{
		using _Base = Base_<UnifiedAllocator<_Ty>, OnGlobal>;
		using _Base::m_rows;
		using _Base::m_cols;
		typedef pointer MyPtr;
	public:
		UnifiedAllocator() : _Base() {}
		template<typename IntegerType>
		UnifiedAllocator(const IntegerType& _M, const IntegerType& _N) :_Base(_M, _N) {}

	private:
		using _Base::m_data;
	};
#pragma endregion
};
}} 
#include "inl\_storage_base.inl"