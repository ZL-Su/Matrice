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
#include <valarray>
#include "_matrix.inl.hpp"
#include "../util/_macros.h"

MATRICE_NAMESPACE_BEGIN_TYPES

/**********************************************************************
						       Matrix view base class 
	    Copyright (c) : Zhilong (Dgelom) Su, since 25/Jul/2018
 **********************************************************************/
template<typename _Ty, typename _Derived> class _View_base 
{
public:
	using value_t = _Ty;
	using value_type = value_t;
	using difference_type = std::ptrdiff_t;
	using pointer = std::add_pointer_t<value_type>;
	using reference = std::add_lvalue_reference_t<value_type>;
	struct range_type
	{
		// _Rang = {from_x, from_y, end_x, end_y} : [from_x, end_x), [from_y, end_y)
		template<typename _Idx, typename = std::enable_if_t<std::is_integral_v<_Idx>>>
		MATRICE_GLOBAL_FINL range_type(const std::initializer_list<_Idx> _Rang)
			:_My_from_x(*_Rang.begin()), _My_from_y(*(_Rang.begin()+1)),
			 _My_end_x(*(_Rang.begin()+2)), _My_end_y(*(_Rang.begin()+3)) {}

		template<typename _Idx, typename = std::enable_if_t<std::is_integral_v<_Idx>>>
		MATRICE_GLOBAL_FINL range_type& operator= (const std::initializer_list<_Idx> _Rang) {
			_My_from_x = *_Rang.begin(), _My_from_y = *(_Rang.begin() + 1);
			_My_end_x = *(_Rang.begin() + 2), _My_end_y = *(_Rang.begin() + 3);
		}

		MATRICE_GLOBAL_FINL auto& begin_x() { return _My_from_x; }
		MATRICE_GLOBAL_FINL auto& begin_y() { return _My_from_y; }
		MATRICE_GLOBAL_FINL auto& end_x() { return _My_end_x; }
		MATRICE_GLOBAL_FINL auto& end_y() { return _My_end_y; }
		MATRICE_GLOBAL_FINL const auto& begin_x() const { return _My_from_x; }
		MATRICE_GLOBAL_FINL const auto& begin_y() const { return _My_from_y; }
		MATRICE_GLOBAL_FINL const auto& end_x() const { return _My_end_x; }
		MATRICE_GLOBAL_FINL const auto& end_y() const { return _My_end_y; }
		MATRICE_GLOBAL_FINL size_t size() const { 
			return(_My_end_x - _My_from_x)*(_My_end_y - _My_from_y); 
		}
		
		difference_type _My_from_x, _My_from_y;
		difference_type _My_end_x, _My_end_y;
	};

	MATRICE_GLOBAL_FINL _View_base(pointer _Ptr, size_t _Size, size_t _Stride, size_t _Offset)
		:_My_data(_Ptr), _My_size(_Size), _My_stride(_Stride), _My_offset(_Offset) {}

	MATRICE_GLOBAL_FINL reference operator[] (size_t i) {
		return _My_data[i*_My_stride];
	}
	MATRICE_GLOBAL_FINL const reference operator[] (size_t i) const {
		return _My_data[i*_My_stride];
	}
	MATRICE_GLOBAL_FINL reference operator() (size_t i) {
		return _My_data[i*_My_stride];
	}
	MATRICE_GLOBAL_FINL const reference operator() (size_t i) const {
		return _My_data[i*_My_stride];
	}

	MATRICE_GLOBAL_FINL auto size() { return (static_cast<_Derived*>(this)->size()); }
	MATRICE_GLOBAL_FINL auto rows() { return (static_cast<_Derived*>(this)->rows()); }
	MATRICE_GLOBAL_FINL auto cols() { return (static_cast<_Derived*>(this)->cols()); }
	MATRICE_GLOBAL_FINL void create(size_t, size_t) {}

	template<typename _Matx>
	MATRICE_GLOBAL_INL auto& operator= (const _Matx& _M) {
		for (size_t i = 0; i < size(); ++i) this->operator()(i) = _M(i);
		return (*static_cast<_Derived*>(this));
	}
	template<typename _Arg> 
	MATRICE_GLOBAL_INL auto& operator= (const Expr::Base_<_Arg>& _Ex) { 
		return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
	template<typename... _Args> 
	MATRICE_GLOBAL_INL auto& operator= (const Expr::MatBinaryExpr<_Args...>& _Ex) { 
		return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
	template<typename... _Args>
	MATRICE_GLOBAL_INL auto& operator= (const Expr::MatUnaryExpr<_Args...>& _Ex) { 
		return (*static_cast<_Derived*>(&_Ex.assign(*this))); }

	template<typename _Op = Expr::EwiseBinaryExpr<_Derived, _Derived, Expr::Op::EwiseSum<value_t>>>
	MATRICE_GLOBAL_FINL auto operator+(const value_t _Right) {
		return _Op(*static_cast<_Derived*>(this), _Right);
	}
	template<typename _Op = Expr::EwiseBinaryExpr<_Derived, _Derived, Expr::Op::EwiseMin<value_t>>>
	MATRICE_GLOBAL_FINL auto operator-(const value_t _Right) {
		return _Op(*static_cast<_Derived*>(this), _Right);
	}
	template<typename _Op = Expr::EwiseBinaryExpr<_Derived, _Derived, Expr::Op::EwiseMul<value_t>>>
	MATRICE_GLOBAL_FINL auto operator*(const value_t _Right) {
		return _Op(*static_cast<_Derived*>(this), _Right);
	}
	template<typename _Op = Expr::EwiseBinaryExpr<_Derived, _Derived, Expr::Op::EwiseDiv<value_t>>>
	MATRICE_GLOBAL_FINL auto operator/(const value_t _Right) {
		return _Op(*static_cast<_Derived*>(this), _Right);
	}
protected:
	pointer _My_data;
	size_t  _My_size;
	size_t  _My_stride;
	size_t  _My_offset;
};

/**********************************************************************
						      Row view for Matrix 
	    Copyright (c) : Zhilong (Dgelom) Su, since 25/Jul/2018
 **********************************************************************/
template<typename _Ty, typename = std::enable_if_t<std::is_arithmetic_v<_Ty>>>
class _Matrix_rview MATRICE_NONHERITABLE : public _View_base<_Ty, _Matrix_rview<_Ty>>
{
	using _Base = _View_base<_Ty, _Matrix_rview<_Ty>>;
	using _Base::_My_data;
	using _Base::_My_size;
	using _Base::_My_stride;
	using _Base::_My_offset;
public:
	using typename _Base::pointer;
	using typename _Base::reference;
	using typename _Base::value_t;
	using _Base::operator+;
	using _Base::operator-;
	using _Base::operator*;
	using _Base::operator/; 
	using _Base::operator=;

	MATRICE_GLOBAL_FINL _Matrix_rview(pointer _Ptr, size_t _Size, size_t _Stride = 1, size_t _Offset = 1)
		:_Base(_Ptr, _Size, _Stride, _Offset) {}

	MATRICE_GLOBAL_FINL auto rows() const { return 1; }
	MATRICE_GLOBAL_FINL auto cols() const { return _My_size; }
	MATRICE_GLOBAL_FINL auto size() const { return _My_size; }

	MATRICE_GLOBAL_FINL value_t sum() const;
};

/**********************************************************************
						     Column view for Matrix 
	    Copyright (c) : Zhilong (Dgelom) Su, since 25/Jul/2018
 **********************************************************************/
template<typename _Ty, typename = std::enable_if_t<std::is_arithmetic_v<_Ty>>>
class _Matrix_cview MATRICE_NONHERITABLE : public _View_base<_Ty, _Matrix_cview<_Ty>>
{
	using _Base = _View_base<_Ty, _Matrix_cview<_Ty>>;
	using _Base::_My_data;
	using _Base::_My_size;
	using _Base::_My_stride;
	using _Base::_My_offset;
public:
	using typename _Base::pointer;
	using typename _Base::reference;
	using typename _Base::value_t;
	using _Base::operator+;
	using _Base::operator-;
	using _Base::operator*;
	using _Base::operator/;
	using _Base::operator=;

	MATRICE_GLOBAL_FINL _Matrix_cview(pointer _Ptr, size_t _Size, size_t _Stride = 1, size_t _Offset = 1)
		:_Base(_Ptr, _Size, _Stride, _Offset) {}

	MATRICE_GLOBAL_FINL auto rows() const { return _My_size; }
	MATRICE_GLOBAL_FINL auto cols() const { return 1; }
	MATRICE_GLOBAL_FINL auto size() const { return _My_size; }

	MATRICE_GLOBAL_FINL value_t sum() const;
};

/**********************************************************************
						      Block view for Matrix 
	    Copyright (c) : Zhilong (Dgelom) Su, since 25/Jul/2018
 **********************************************************************/
template<typename _Ty, typename = std::enable_if_t<std::is_arithmetic_v<_Ty>>>
class _Matrix_block MATRICE_NONHERITABLE : public _View_base<_Ty, _Matrix_block<_Ty>>
{
	using _Base = _View_base<_Ty, _Matrix_block<_Ty>>;
	using typename _Base::range_type;
	using typename _Base::difference_type;
	using _Base::_My_data;   //begin of this block data
	using _Base::_My_size;   //cols of this block
	using _Base::_My_stride; //cols of source matrix
	using _Base::_My_offset; //offset relative to original matrix data
public:
	using typename _Base::pointer;
	using typename _Base::reference;
	using typename _Base::value_t;
	using _Base::operator+;
	using _Base::operator-;
	using _Base::operator*;
	using _Base::operator/;
	using _Base::operator=;

	MATRICE_GLOBAL_FINL _Matrix_block(pointer _Ptr, size_t _Cols, const range_type _Range)
		: _Base(_Ptr + _Range.begin_x() + _Range.begin_y()*_Cols, 
			_Range.end_x() - _Range.begin_x(), _Cols, 
			_Range.begin_x() + _Range.begin_y()*_Cols),
		_My_range(_Range) {}

	//i zero-based local row index
	MATRICE_GLOBAL_FINL pointer operator[] (difference_type i) {
		return (_My_data + i * _My_stride);
	}
	MATRICE_GLOBAL_FINL const pointer operator[] (difference_type i) const {
		return (_My_data + i * _My_stride);
	}
	//i zero-based local linear index
	MATRICE_GLOBAL_FINL reference operator() (difference_type i) {
		auto _Row = i / size();
		return this->operator[](_Row)[i - _My_size * _Row];
	}
	MATRICE_GLOBAL_FINL const reference operator() (difference_type i) const {
		auto _Row = i / size();
		return this->operator[](_Row)[i - _My_size * _Row];
	}

	MATRICE_GLOBAL_FINL auto rows() const { return _My_range.end_y() - _My_range.begin_y(); }
	MATRICE_GLOBAL_FINL auto cols() const { return _My_size; }
	MATRICE_GLOBAL_FINL auto size() const { return rows()*cols(); }

private:
	range_type _My_range;
};
MATRICE_NAMESPACE_END_TYPES
