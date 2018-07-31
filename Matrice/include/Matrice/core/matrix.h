/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once
#include <type_traits>
#include <array>
#include "../private/_matrix_base.hpp"
#include "../private/_dev_matrix_base.h"

MATRICE_NAMESPACE_BEGIN_TYPES

/*******************************************************************
	Generic Matrix Class with Aligned Static Memory Allocation
	          Copyright (c) : Zhilong Su 14/Feb/2018
*******************************************************************/
template<typename _Ty, int _M, int _N>
class Matrix_ : public Base_<_Ty, _M, _N>
{
	typedef Base_<_Ty, _M, _N> base_t;
	typedef Matrix_                                     Myt;
	typedef const Myt&                         const_my_ref;
	using base_t::m_rows;
	using base_t::m_cols;
	using base_t::m_data;
public:
	using typename base_t::value_t;
	using typename base_t::pointer;
	using typename base_t::const_init_list;
	enum { Size = _M * _N, CompileTimeRows = _M, CompileTimeCols = _N, };
	MATRICE_GLOBAL_FINL Matrix_(const value_t _val) noexcept : base_t({_val}) {};
	MATRICE_GLOBAL_FINL Matrix_(pointer data) noexcept : base_t(_M, _N, data) {};
	MATRICE_HOST_FINL   Matrix_(const_init_list _list) noexcept : base_t(_list) {}
	MATRICE_GLOBAL_FINL Matrix_(Myt&& _other) noexcept : base_t(_other) {}
	MATRICE_GLOBAL_FINL Matrix_(const_my_ref _other) noexcept : base_t(_other) {};
	MATRICE_GLOBAL_FINL Matrix_(int _pld1=0, int _pld2=0) noexcept : base_t() {};
	template<typename _Arg> MATRICE_GLOBAL_FINL Matrix_(const _Arg& _arg) noexcept : base_t(_arg) {};

	MATRICE_HOST_FINL Myt& operator= (const_init_list _list) { return base_t::operator=(_list); }
	MATRICE_GLOBAL_FINL Myt& operator= (const_my_ref _other) { return base_t::operator=(_other); }
	MATRICE_GLOBAL_FINL Myt& operator= (Myt&& _other) { return base_t::operator=(std::move(_other)); }
	template<typename _Arg> MATRICE_GLOBAL_FINL Myt& operator= (const _Arg& _arg) { return base_t::operator=(_arg); }
	MATRICE_GLOBAL_FINL operator std::array<value_t, Size>() const { return _Fill_array<value_t, Size>(base_t::begin()); }
	MATRICE_GLOBAL_FINL operator std::array<value_t, Size>() { return _Fill_array<value_t, Size>(base_t::begin()); }
	MATRICE_GLOBAL_FINL operator Matrix_<value_t, __, __>() const { return Matrix_<value_t, __, __>(m_rows, m_cols, m_data); }
	MATRICE_GLOBAL_FINL operator Matrix_<value_t, __, __>() { return Matrix_<value_t, __, __>(m_rows, m_cols, m_data); }
#ifdef __use_ocv_as_view__
	using base_t::operator ocv_view_t;
#endif
};
//matrix with host managed memory allocator [float]
template<int _Rows, int _Cols> using Matrixf = Matrix_<float, 
	compile_time_size<_Rows, _Cols>::CompileTimeRows, 
	compile_time_size<_Rows, _Cols>::CompileTimeCols>;
//matrix with host managed memory allocator [double]
template<int _Rows, int _Cols> using Matrixd = Matrix_<double,
	compile_time_size<_Rows, _Cols>::CompileTimeRows,
	compile_time_size<_Rows, _Cols>::CompileTimeCols>;

/*******************************************************************
    Generic Matrix Class with Aligned Dynamic Memory Allocation
	         Copyright (c) : Zhilong Su 14/Feb/2018
 ******************************************************************/
template<typename _Ty>
class Matrix_<_Ty, __, __> : public Base_<_Ty, __, __>
{
	typedef Base_<_Ty, __, __> base_t;
	typedef Matrix_                                 Myt;
	typedef const Myt&                     const_my_ref;
	using base_t::m_data;
	using base_t::m_rows;
	using base_t::m_cols;
public:
	using typename base_t::value_t;
	using typename base_t::pointer;
	using typename base_t::const_init_list;
	enum { Size = __, CompileTimeRows = __, CompileTimeCols = __, };
	MATRICE_GLOBAL_FINL Matrix_(int _rows) noexcept : base_t(_rows, 1) {};
	MATRICE_GLOBAL_FINL Matrix_(Myt&& _other) noexcept : base_t(_other) {};
	MATRICE_GLOBAL_FINL Matrix_(const_my_ref _other) noexcept : base_t(_other) {};
	template<typename... _Args> MATRICE_GLOBAL_FINL Matrix_(const _Args&... args) noexcept : base_t(args...) {};
	//MATRICE_GLOBAL_INL Matrix_() noexcept : base_t() {};
	//MATRICE_GLOBAL_INL Matrix_(int _rows, int _cols) noexcept : base_t(_rows, _cols) {};
	//MATRICE_GLOBAL_INL Matrix_(int _rows, int _cols, pointer _data) noexcept : base_t(_rows, _cols, _data) {};
	//MATRICE_GLOBAL_INL Matrix_(int _rows, int _cols, const value_t _val) noexcept : base_t(_rows, _cols, _val) {};
	//template<int _M, int _N>
	//MATRICE_GLOBAL_INL Matrix_(const Matrix_<value_t, _M, _N>& _other) noexcept : base_t(_other.rows(), _other.cols(), _other.data()) {};
	//template<typename _Expr> MATRICE_GLOBAL_INL Matrix_(const _Expr& _other) noexcept : base_t(_other) {};

	template<typename _Arg> MATRICE_GLOBAL_FINL Myt& operator= (const _Arg& _arg) { return base_t::operator=(_arg); }
	MATRICE_GLOBAL_FINL Myt& operator= (Myt&& _other) { return base_t::operator= (std::move(_other)); }
	MATRICE_GLOBAL_FINL Myt& operator= (const_my_ref& _other) { return base_t::operator=(_other); }
	MATRICE_HOST_FINL Myt& operator= (const_init_list _list) { return base_t::operator=(_list); }
	MATRICE_GLOBAL void create(int_t rows, int_t cols = 1);
	MATRICE_GLOBAL void create(int_t rows, int_t cols, value_t _val);
};
//matrix with host dynamic memory allocator
template<typename _Ty> using Matrix = Matrix_<_Ty,
	compile_time_size<>::RunTimeDeducedInHost,
	compile_time_size<>::RunTimeDeducedInHost>;

#if (defined __enable_cuda__ && !defined __disable_cuda__)
/*******************************************************************
    Generic Matrix Class with Unifief Memory Allocation (CUDA)
	         Copyright (c) : Zhilong Su 24/May/2018
 ******************************************************************/
template<typename _Ty>
class Matrix_<_Ty, -1, __> : public Base_<_Ty, -1, __>
{
	typedef Base_<_Ty, -1, __>                   base_t;
	typedef Matrix_                                 Myt;
	typedef const Myt&                     const_my_ref;
	using base_t::m_data;
	using base_t::m_rows;
	using base_t::m_cols;
public:
	using typename base_t::value_t;
	using typename base_t::pointer;
	using typename base_t::const_init_list;
	enum { Size = __, CompileTimeRows = -1, CompileTimeCols = __, };
	MATRICE_GLOBAL_INL Matrix_(int _rows) noexcept : base_t(_rows, 1) {};
	MATRICE_GLOBAL_INL Matrix_(Myt&& _other) noexcept : base_t(_other) {};
	template<typename... _Args> MATRICE_GLOBAL_INL Matrix_(_Args... args) noexcept : base_t(args...) {};

	template<typename _Arg> MATRICE_GLOBAL_INL Myt& operator= (const _Arg& _arg) { return base_t::operator=(_arg); }
	MATRICE_GLOBAL_INL Myt& operator= (Myt&& _other) { return base_t::operator=(std::move(_other)); }
	MATRICE_HOST_INL Myt& operator= (const_init_list _list) { return base_t::operator=(_list); }
	MATRICE_GLOBAL void create(int_t rows, int_t cols = 1);
};
//matrix with CUDA unified memory allocator
template<typename _Ty> using Umatrix = Matrix_<_Ty,
	compile_time_size<>::RunTimeDeducedInDevice,
	compile_time_size<>::RunTimeDeducedInHost>;

/*******************************************************************
    Generic Matrix Class with Device Memory Allocation (CUDA)
	         Copyright (c) : Zhilong Su 25/May/2018
 ******************************************************************/
template<typename _Ty>
class Matrix_<_Ty, -1, -1> : public Base_<_Ty, -1, -1>, public device::Base_<_Ty>
{
	typedef Base_<_Ty, -1, -1>                   base_t;
	typedef device::Base_<_Ty>            device_base_t;
	typedef Matrix_                                 Myt;
	typedef const Myt&                     const_my_ref;
	using base_t::m_data;
	using base_t::m_rows;
	using base_t::m_cols;
	std::size_t m_pitch = base_t::m_storage.pitch();
public:
	enum { Size = -1, CompileTimeRows = -1, CompileTimeCols = -1, };
	using typename base_t::value_t;
	using typename base_t::pointer;
	using typename base_t::const_init_list;
	using device_base_t::operator+;

	MATRICE_GLOBAL_INL Matrix_(int _rows) noexcept : base_t(_rows, 1), 
		device_base_t(m_data, &m_pitch, &m_cols, &m_rows) {};
	MATRICE_GLOBAL_INL Matrix_(const Myt& _other) noexcept : base_t(_other),
		device_base_t(m_data, &m_pitch, &m_cols, &m_rows) {};
	MATRICE_GLOBAL_INL Matrix_(Myt&& _other) noexcept : base_t(_other), 
		device_base_t(m_data, &m_pitch, &m_cols, &m_rows) {};
	template<int _M = 0, int _N = _M>
	MATRICE_HOST_INL Matrix_(const Matrix_<_Ty, _M, _N>& _other) noexcept: base_t(_other),
		device_base_t(m_data, &m_pitch, &m_cols, &m_rows) {};
	template<typename... _Args> MATRICE_GLOBAL_INL Matrix_(const _Args&... args) noexcept : base_t(args...), 
		device_base_t(m_data, &m_pitch, &m_cols, &m_rows) {};

	template<typename _Arg> 
	MATRICE_GLOBAL_INL Myt& operator= (const _Arg& _arg) { return device_base_t::operator=(_arg); }
	MATRICE_DEVICE_INL Myt& operator= (const Myt& _other) { return device_base_t::operator=(_other); }
	MATRICE_DEVICE_INL Myt& operator= (Myt&& _other) { return device_base_t::operator=(std::move(_other)); }
	MATRICE_GLOBAL_INL Myt& operator= (const_init_list _list) { return device_base_t::operator=(_list.begin()); }
	/*template<int _M, int _N>
	MATRICE_GLOBAL_INL operator Matrix_<value_t, _M, _N>() const {
		Matrix_<value_t, _M, _N> _Ret(m_rows, m_cols);
		privt::unified_sync<value_t, base_t::m_storage.location, _Ret.location, base_t::m_storage.option>::op(_Ret.data(), m_data, m_rows, m_cols, m_pitch);
		return std::move(_Ret);
	}*/
	MATRICE_GLOBAL void create(int_t rows, int_t cols = 1);

};
//matrix with CUDA device memory allocator
template<typename _Ty> using Dmatrix = Matrix_<_Ty,
	compile_time_size<>::RunTimeDeducedInDevice,
	compile_time_size<>::RunTimeDeducedInDevice>;
#endif

#pragma region Matrix base which using storage class as template arg. 
template<typename _Ty, class _Storage = details::Storage_<_Ty>::
#ifdef __CXX11__
	SharedAllocator<_Ty>
#else
	DynamicAllocator<_Ty>
#endif
>
class MatrixBase_
#ifdef __use_ocv_as_view__
	: public ocv_view_t
#endif
{
	typedef typename details::Storage_<_Ty>::value_t   value_t;
	typedef typename details::Storage_<_Ty>::pointer   pointer;
	typedef typename details::Storage_<_Ty>::reference reference;
	typedef typename details::Storage_<_Ty>::idx_t     idx_t;
	typedef typename Location                          loctn_t;
	typedef MatrixBase_ Myt;
public:
	explicit MatrixBase_(const int_t _m, const int_t _n) noexcept;
	explicit MatrixBase_() noexcept 
		: m_storage(), m_rows(m_storage.rows()), m_cols(m_storage.cols())
	{
		if (m_storage.location() == loctn_t::OnStack)
		{
			m_data = m_storage.data();
#ifdef __use_ocv_as_view__
			*static_cast<ocv_view_t*>(this) = ocv_view_t(m_storage.rows(), m_storage.cols(), ocv_view_t_cast<value_t>::type, m_data);
#endif
		}
	};

public:
	inline value_t& operator[] (int_t idx) const { return m_data[idx]; }
	inline value_t& operator() (int_t ridx, int_t cidx) const { return m_data[cidx+ridx* m_cols]; }
	inline pointer ptr(int_t ridx) const { return (m_data + ridx*m_cols); }
	inline value_t& colmaj(int_t ridx, int_t cidx) const { return m_data[ridx + cidx * m_rows]; }

protected:
	pointer m_data;
	_Storage m_storage;
private:
	std::size_t m_rows, m_cols;
};
#pragma endregion
MATRICE_NAMESPACE_END_TYPES
namespace dge = dgelom::types;