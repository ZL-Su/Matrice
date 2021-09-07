/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
#include <iterator>
#include <vector>
#include "util/_macros.h"
#include "_shape.hpp"

MATRICE_NAMESPACE_BEGIN_

template<typename _InIt> MATRICE_GLOBAL_FINL
_InIt _End(const _InIt _Begin, size_t _Size, size_t _Stride = 1) {
	return (_Begin + _Size*_Stride);
}

/**********************************************************************
    Forward iterator for range [_Myptr, _Myend), which is compatible 
                         with STD::ITERATOR
	    Copyright (c) : Zhilong (Dgelom) Su, since 12/Jul/2018
 **********************************************************************/
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
class _Iterator {
	using _Myt = _Iterator;
public:
	using iterator_category = std::random_access_iterator_tag;
	using value_type = _Ty;
	using pointer = value_type*;
	using reference = value_type&;
	using difference_type = std::ptrdiff_t;
	enum { rows_at_compiletime = 0, cols_at_compiletime = 0 };

	MATRICE_GLOBAL_FINL _Iterator(pointer _Ptr) noexcept 
		:_Myptr(_Ptr), _Mybegin(_Ptr), _Mysize(0), _Mystep(1) {}
	MATRICE_GLOBAL_FINL _Iterator(pointer _Ptr, size_t _Size, size_t _Step =  1) noexcept 
		:_Myptr(_Ptr), _Mybegin(_Ptr), _Mysize(_Size), _Mystep(_Step) {}
	
	MATRICE_GLOBAL_FINL const reference operator*() const noexcept { 
		return ((reference)(*_Myptr));
	}
	MATRICE_GLOBAL_FINL reference operator*() noexcept {
		return ((reference)(*_Myptr));
	}
	MATRICE_GLOBAL_FINL pointer operator->() const { 
		return (std::pointer_traits<pointer>::pointer_to(**this)); 
	}
	MATRICE_GLOBAL_FINL _Myt& operator++() { //preincrement
		_Myptr += _Mystep;
		return (*this);
	}
	MATRICE_GLOBAL_FINL _Myt operator++(int) noexcept { //postincrement
		auto _Tmp = *this;
		*this += _Mystep;
		return (_Tmp);
	}
	MATRICE_GLOBAL_FINL _Myt& operator--() noexcept { //preincrement
		_Myptr -= _Mystep;
		return (*this);
	}
	MATRICE_GLOBAL_FINL _Myt operator--(int) noexcept { //postincrement
		auto _Tmp = *this;
		*this += _Mystep;
		return (_Tmp);
	}
	MATRICE_GLOBAL_FINL _Myt& operator+=(difference_type _Offset) {
		_Offset *= _Mystep;
#if _ITERATOR_DEBUG_LEVEL == 2
		if (_Offset != 0) {
			if (_Myptr + _Offset < _Myptr || _Myend < _Myptr + _Offset) {
				DGELOM_ERROR("iterator + offset out of range");
			}
		}
#endif
		_Myptr += _Offset;
		return (*this);
	}
	MATRICE_GLOBAL_FINL _Myt operator+(difference_type _Offset) const {
		auto _Tmp = *this;
		return (_Tmp += _Offset);
	}
	MATRICE_GLOBAL_FINL _Myt& operator-=(difference_type _Offset) {
		return (*this += -(_Offset * _Mystep));
	}
	MATRICE_GLOBAL_FINL _Myt operator-(difference_type _Offset) const noexcept {
		auto _Tmp = *this;
		return (_Tmp -= (_Offset));
	}
	MATRICE_GLOBAL_FINL difference_type operator-(const _Myt& _Right) const noexcept {
		return (_Myptr - _Right._Myptr);
	}
	MATRICE_GLOBAL_FINL reference operator[](difference_type _Offset) const noexcept {
		return (*(*this + _Offset * _Mystep));
	}
	MATRICE_GLOBAL_FINL bool operator==(const _Myt& _Right) const noexcept {
		return (_Myptr == _Right._Myptr);
	}
	MATRICE_GLOBAL_FINL bool operator!=(const _Myt& _Right) const noexcept {
		return (!(*this == _Right));
	}
	MATRICE_GLOBAL_FINL bool operator<(const _Myt& _Right) const noexcept {
		return (_Myptr < _Right._Myptr);
	}
	MATRICE_GLOBAL_FINL bool operator>(const _Myt& _Right) const noexcept {
		return (_Right < *this);
	}
	MATRICE_GLOBAL_FINL bool operator<=(const _Myt& _Right) const noexcept {
		return (!(_Right < *this));
	}
	MATRICE_GLOBAL_FINL bool operator>=(const _Myt& _Right) const noexcept {
		return (!(_Right > *this));
	}

	// \test for iterator end condition
	MATRICE_GLOBAL_FINL operator bool() const noexcept { return (_Myptr != _Myend); }

	// \return pointer to current object
	MATRICE_GLOBAL_FINL operator pointer() noexcept { return (_Myptr); }

	// \forward range iteration methods for [this->_Myptr, this->_Myend)
	MATRICE_GLOBAL_FINL decltype(auto) begin() noexcept {
		return (*this); 
	}
	MATRICE_GLOBAL_FINL auto end() noexcept {
		auto _Tmp = *this; 
		return (_Tmp += _Mysize); 
	}
	MATRICE_GLOBAL_FINL decltype(auto) begin() const noexcept {
		return (*this); 
	}
	MATRICE_GLOBAL_FINL const auto end() const noexcept {
		auto _Tmp = *this; 
		return (_Tmp += _Mysize);
	}

	/// <summary>
	/// \brief Get distance between adjacent iterator positions.
	/// </summary>
	MATRICE_GLOBAL_FINL decltype(auto)stride()const noexcept {
		return (_Mystep);
	}

	// \current iterator position
	MATRICE_GLOBAL_FINL size_t pos() const noexcept { 
		return std::distance(_Mybegin, _Myptr)/_Mystep; 
	}

protected:
	size_t _Mysize;
	size_t _Mystep;
	pointer _Myptr = nullptr;
	pointer _Myend = _End(_Myptr, _Mysize, _Mystep);
	pointer _Mylast = _Myend - _Mystep;

private:
	const pointer _Mybegin = nullptr;
};
template<typename _Ty>
MATRICE_GLOBAL_FINL _Iterator<_Ty> operator+ (typename _Iterator<_Ty>::difference_type _Offset, _Iterator<_Ty> _Next) {
	return (_Next += _Offset);
}

/**********************************************************************
						Matrix Const Iterator
		Copyright (c) : Zhilong (Dgelom) Su, since May/23/2020
 **********************************************************************/
template<typename _Ty>
class _Matrix_const_iterator : public _Iterator<_Ty>
{
	using _Mybase = _Iterator<_Ty>;
	using _Myt = _Matrix_const_iterator;
public:
	using typename _Mybase::pointer;
	using typename _Mybase::value_type;
	using typename _Mybase::reference;
	_Matrix_const_iterator(pointer ptr=nullptr, size_t size=0) noexcept 
		:_Mybase{ ptr, size } {
	}

	MATRICE_GLOBAL_FINL _Myt& operator=(const _Myt& other) noexcept {
		_Mybase::_Mysize = other._Mysize;
		_Mybase::_Myptr = other._Myptr;
		_Mybase::_Myend = other._Myend;
		_Mybase::_Mystep = other._Mystep;
		_Mybase::_Mylast = other._Mylast;
		return (*this);
	}
	
	MATRICE_GLOBAL_FINL const reference operator*() const noexcept {
		return _Mybase::operator*();
	}
	MATRICE_GLOBAL_FINL _Myt& operator++() noexcept {
		_Mybase::_Myptr += _Mybase::_Mystep;
		return (*this);
	}
	MATRICE_GLOBAL_FINL _Myt& operator++(int) noexcept {
		auto _Pre = (*this);
		++(*this);
		return (_Pre);
	}
	MATRICE_GLOBAL_FINL _Myt& operator--() noexcept {
		_Mybase::_Myptr -= _Mybase::_Mystep;
		return (*this);
	}
	MATRICE_GLOBAL_FINL _Myt& operator--(int) noexcept {
		auto _Pre = (*this);
		--(*this);
		return (_Pre);
	}
	MATRICE_GLOBAL_FINL operator pointer() noexcept {
		return _Mybase::_Myptr;
	}
	MATRICE_GLOBAL_FINL _Myt& stride(size_t step) noexcept {
		_Mybase::_Mystep = step;
		_Mybase::_Myend = _End(_Myptr, _Mybase::_Mysize, _Mybase::_Mystep);
		return (*this);
	}
	using _Mybase::operator+=;
	using _Mybase::operator-=;
	using _Mybase::operator!=;
	using _Mybase::operator==;
	using _Mybase::operator<;
	using _Mybase::operator<=;
	using _Mybase::operator>;
	using _Mybase::operator>=;
	using _Mybase::operator[];
	using _Mybase::operator+;
	using _Mybase::operator-;
	using _Mybase::operator bool;
	using _Mybase::pos;

private:
	using _Mybase::_Myptr;
};

/**********************************************************************
						Forward Range Iterator
	    Copyright (c) : Zhilong (Dgelom) Su, since 12/Jul/2018
 **********************************************************************/
template<typename _Ty>
class _Matrix_forward_iterator : public _Iterator<_Ty>
{
	using _Mybase = _Iterator<_Ty>;
public:
	template<typename... _Args> 
	MATRICE_GLOBAL_FINL _Matrix_forward_iterator(_Args... _args)noexcept
		: _Mybase(_args...) {}

private:
	typename _Mybase::difference_type _Myidx, _Myidy;
};

/**********************************************************************
					Row-wise Forward Range Iterator
	    Copyright (c) : Zhilong (Dgelom) Su, since 12/Jul/2018
 **********************************************************************/
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
class _Matrix_rwise_iterator : public _Iterator<_Ty>
{
	/*<Note> 
		_Myptr  := input pointer to row
		_Mysize := number of rows
		_Mystep := number of cols
		_Return  := i-th element in the row that _Myptr point to
	  </Note>*/
	using _Mybase = _Iterator<_Ty>;
	using _Myt = _Matrix_rwise_iterator;
public:
	template<typename... _Args>
	MATRICE_GLOBAL_FINL _Matrix_rwise_iterator(_Args... _args) noexcept
		: _Mybase(_args...) {}

	// \element-wise accessor
	MATRICE_GLOBAL_FINL auto& operator[](typename _Mybase::difference_type _Offset) const {
		return (*(_Mybase::_Myptr + _Offset));
	}
	// \copy data from another iterator
	MATRICE_GLOBAL_FINL _Myt& operator= (const _Mybase& _Iter) {
		std::copy(_Iter.begin(), _Iter.end(), _Mybase::_Myptr);
		return (*this);
	}
	// \copy data from another pointer
	MATRICE_GLOBAL_FINL _Myt& operator= (const typename _Mybase::pointer _Ptr) {
		std::copy(_Ptr, _Ptr + _Mybase::_Mystep, _Mybase::_Myptr);
		return (*this);
	}
	// \copy data from initializer list
	MATRICE_GLOBAL_FINL _Myt& operator= (const std::initializer_list<_Ty> _List) {
		std::copy(_List.begin(), _List.end(), _Mybase::_Myptr);
		return (*this);
	}
	// \iterator to tranverse the element in current row 
	MATRICE_GLOBAL_FINL auto begin() {
		return _Mybase(_Mybase::_Myptr, _Myrange);
	}
	MATRICE_GLOBAL_FINL const auto begin() const {
		return _Mybase(_Mybase::_Myptr, _Myrange);
	}
	MATRICE_GLOBAL_FINL auto end() {
		return _Mybase(_Mybase::_Myptr + _Myrange, _Myrange);
	}
	MATRICE_GLOBAL_FINL const auto end() const {
		return _Mybase(_Mybase::_Myptr + _Myrange, _Myrange);
	}
private:
	size_t _Myrange = _Mybase::_Mystep;
};

/**********************************************************************
					 Column-wise Forward Range Iterator
	    Copyright (c) : Zhilong (Dgelom) Su, since 12/Jul/2018
 **********************************************************************/
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
class _Matrix_cwise_iterator : public _Iterator<_Ty>
{
	/*<Note> 
		_Myptr  := input pointer to column
		_Mysize := number of cols
		_Mystep := 1
		_Return := i-th element in the column that _Myptr point to
	  </Note>*/
	using _Mybase = _Iterator<_Ty>;
	using _Myt = _Matrix_cwise_iterator;
public:
	template<typename... _Args>
	MATRICE_GLOBAL_FINL _Matrix_cwise_iterator(_Args... _args) noexcept
		: _Mybase(_args...) {
		std::swap(_Myrange, _Mybase::_Mystep);
	}
	// \element-wise accessor
	MATRICE_GLOBAL_FINL auto& operator[](typename _Mybase::difference_type _Offset) const {
		return (*(_Mybase::_Myptr + _Offset*_Mybase::_Mysize));
	}
	// \copy data from another iterator
	MATRICE_GLOBAL_FINL _Myt& operator= (const _Mybase& _Iter) {
		_Mybase _Tmp(_Mybase::_Myptr, _Myrange, _Mybase::_Mysize);
		std::copy(_Iter.begin(), _Iter.end(), _Tmp.begin());
		return (*this);
	}
	// \copy data from another pointer
	MATRICE_GLOBAL_FINL _Myt& operator= (const typename _Mybase::pointer _Ptr) {
		_Mybase _Tmp(_Mybase::_Myptr, _Myrange, _Mybase::_Mysize);
		std::copy(_Ptr, _Ptr + _Mybase::_Mystep, _Tmp.begin());
		return (*this);
	}
	// \copy data from initializer list
	MATRICE_GLOBAL_FINL _Myt& operator= (const std::initializer_list<_Ty> _List) {
		_Mybase _Tmp(_Mybase::_Myptr, _Myrange, _Mybase::_Mysize);
		std::copy(_List.begin(), _List.end(), _Tmp.begin());
		return (*this);
	}
	// \iterator to tranverse the element in current column
	MATRICE_GLOBAL_FINL auto begin() {
		return _Mybase(_Mybase::_Myptr, _Myrange, _Mybase::_Mysize);
	}
	MATRICE_GLOBAL_FINL const auto begin() const {
		return _Mybase(_Mybase::_Myptr, _Myrange, _Mybase::_Mysize);
	}
	MATRICE_GLOBAL_FINL auto end() {
		return _Mybase(_End(_Mybase::_Myptr, _Myrange, _Mybase::_Mysize), _Myrange, _Mybase::_Mysize);
	}
	MATRICE_GLOBAL_FINL const auto end() const {
		return _Mybase(_End(_Mybase::_Myptr, _Myrange, _Mybase::_Mysize), _Myrange, _Mybase::_Mysize);
	}
private:
	size_t _Myrange = 1;
};


/**********************************************************************
			    Forward Grid Iterator[Not Yet Completed]
		Copyright (c) : Zhilong (Dgelom) Su, since Apr/26/2020
 **********************************************************************/
template<typename _Ty>
class _Matrix_grid_iterator
{
	using _Myt = _Matrix_grid_iterator;
public:
	using iterator_category = std::random_access_iterator_tag;
	using value_type = _Ty;
	using pointer = std::add_pointer_t<value_type>;
	using reference = std::add_lvalue_reference_t<typename std::pointer_traits<pointer>::element_type>;
	using difference_type = std::ptrdiff_t;
	enum { rows_at_compiletime = 0, cols_at_compiletime = 0 };

	MATRICE_GLOBAL_FINL _Matrix_grid_iterator(pointer _Ptr, shape_t<3> _Shape, diff_t _Incy = 1, diff_t _Incx = 1) noexcept
		:_Mybegin(_Ptr), _Myptr(_Ptr), _Myshape(_Shape),
		_Myincx(_Incx), _Myincy(_Incy), _Mysize(_Shape.size()), 
		_Myend(_Ptr+_Shape.size()){
	}

	MATRICE_GLOBAL_FINL reference operator*() const {
		return (reference(*_Myptr));
	}

	MATRICE_GLOBAL_FINL _Myt& operator++() { //preincrement
		_Mycol += _Myincx;
		if (_Mycol >= _Myshape.cols()) {
			_Mycol %= _Myshape.cols();
			_Myrow += _Myincy;
		}
		_Myptr = _Mybegin + _Mycol + _Myrow * _Myshape.cols();
		return (*this);
	}

	MATRICE_GLOBAL_FINL bool operator==(const _Myt& _Right) const {
		return (_Myptr == _Right._Myptr);
	}
	MATRICE_GLOBAL_FINL bool operator!=(const _Myt& _Right) const {
		return (!(*this == _Right));
	}

	MATRICE_GLOBAL_FINL decltype(auto) begin() noexcept { 
		return (*this); 
	}
	MATRICE_GLOBAL_FINL decltype(auto) begin() const noexcept { 
		return (*this); 
	}
	MATRICE_GLOBAL_FINL const auto end() const noexcept {
		auto _Tmp = *this;
		//auto _Col = _Mycol + _Myincx, _Row = _Myrow;
		//if (_Col >= _Myshape.cols()) {
		//	_Col = 0;//%= _Myshape.cols();
		//	_Row += _Myincy;
		//}
		return (++_Tmp);
	}

	MATRICE_GLOBAL_FINL operator bool() const noexcept {
		return _Myptr < _Myend;
	}

	/**
	 *\brief Get the position of the current element the iterator points to.
	 *\ex const auto[row, cols] = it.pos();
	 */
	MATRICE_GLOBAL_FINL auto pos() const noexcept {
		return std::make_tuple(_Myrow, _Mycol);
	}

	MATRICE_GLOBAL_FINL _Myt& reset_strides(diff_t _Incx, diff_t _Incy)noexcept {
		_Myincx = _Incx, _Myincy = _Incy;
		return (*this);
	}
private:
	const pointer _Mybegin = nullptr, _Myend = nullptr;
	pointer _Myptr = nullptr;
	shape_t<3> _Myshape;
	size_t _Mycol = 0, _Myrow = 0, _Mysize = 0;
	difference_type _Myincx, _Myincy;
};

/**********************************************************************
						 Transform Forward Iterator
		 Copyright (c) : Zhilong (Dgelom) Su, since May/10/2019
 **********************************************************************/
template<typename _It, typename _Op>
class _Transform_iterator {
	using _Myt = _Transform_iterator;
public:
	using iterator_category = conditional_t<is_pointer_v<_It>, std::random_access_iterator_tag, typename _It::iterator_category>;
	using value_type = decltype(std::declval<_Op>()(*std::declval<_It>()));
	using pointer = std::add_pointer_t<value_type>;
	using reference = std::add_lvalue_reference_t<value_type>;
	using difference_type = std::ptrdiff_t;
	enum { rows_at_compiletime = 0, cols_at_compiletime = 0 };

	MATRICE_HOST_INL _Transform_iterator(const _It& it, _Op&& op)
		: _Myit(it), _Myop(op) {}
	MATRICE_HOST_INL _Transform_iterator(const _Myt& other)
		: _Myit(other._Myit), _Myop(other._Myop) {}

	MATRICE_HOST_INL _Myt& operator++() { ++_Myit; return (*this); }
	MATRICE_HOST_INL difference_type operator-(const _Myt& other) const { 
		return(_Myit - other._Myit); 
	}
	MATRICE_HOST_INL value_type operator*() const { return (_Myop(*_Myit)); }
	MATRICE_HOST_INL bool operator!=(const _Myt& other) const { 
		return(_Myit != other._Myit); 
	}

private:
	_It _Myit;
	_Op _Myop;
};
template<typename _It, typename _Op>
MATRICE_HOST_INL _Transform_iterator<_It,_Op> make_transform_iter(const _It& it, _Op&& op) {
	return _Transform_iterator<_It, _Op>(it, op);
}
_MATRICE_NAMESPACE_END
