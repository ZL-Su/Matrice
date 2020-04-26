/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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
#include "../util/_macros.h"

MATRICE_NAMESPACE_BEGIN_

template<typename _InIt> MATRICE_GLOBAL_FINL
_InIt _End(const _InIt _Begin, size_t _Size, size_t _Stride = 1) {
	return (_Begin + _Size*_Stride);
}

/**********************************************************************
 Base for Forward Range Iterator[_Myptr, _Myend), which is compatible 
                         with STD::ITERATOR
	    Copyright (c) : Zhilong (Dgelom) Su, since 12/Jul/2018
 **********************************************************************/
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
class iterator_base {
public:
	using iterator_category = std::random_access_iterator_tag;
	using value_type = _Ty;
	using pointer = std::add_pointer_t<value_type>;
	using reference = std::add_lvalue_reference_t<typename std::pointer_traits<pointer>::element_type>;
	using difference_type = std::ptrdiff_t;
	enum { rows_at_compiletime = 0, cols_at_compiletime = 0 };

	MATRICE_GLOBAL_FINL iterator_base(pointer _Ptr) noexcept 
		:_Myptr(_Ptr), _Mybegin(_Ptr), _Mysize(0), _Mystep(0) {}
	MATRICE_GLOBAL_FINL iterator_base(pointer _Ptr, size_t _Size, size_t _Step =  1) noexcept 
		:_Myptr(_Ptr), _Mybegin(_Ptr), _Mysize(_Size), _Mystep(_Step) {}
	
	MATRICE_GLOBAL_FINL reference operator*() const { 
		return ((reference)(*_Myptr));
	}
	MATRICE_GLOBAL_FINL pointer operator->() const { 
		return (std::pointer_traits<pointer>::pointer_to(**this)); 
	}
	MATRICE_GLOBAL_FINL iterator_base operator++() { //preincrement
		_Myptr += _Mystep;
		return (*this);
	}
	MATRICE_GLOBAL_FINL iterator_base operator++(int) { //postincrement
		auto _Tmp = *this;
		*this += _Mystep;
		return (_Tmp);
	}
	MATRICE_GLOBAL_FINL iterator_base operator--() { //preincrement
		_Myptr -= _Mystep;
		return (*this);
	}
	MATRICE_GLOBAL_FINL iterator_base operator--(int) { //postincrement
		auto _Tmp = *this;
		*this += _Mystep;
		return (_Tmp);
	}
	MATRICE_GLOBAL_FINL iterator_base& operator+=(difference_type _Offset) {
		_Offset *= _Mystep;
#if _ITERATOR_DEBUG_LEVEL == 2
		if (_Offset != 0) {
			if (_Myptr + _Offset < _Myptr || _Myend < _Myptr + _Offset) {
				throw std::runtime_error("iterator + offset out of range");
			}
		}
#endif
		_Myptr += _Offset;
		return (*this);
	}
	MATRICE_GLOBAL_FINL iterator_base operator+(difference_type _Offset) const {
		auto _Tmp = *this;
		return (_Tmp += _Offset);
	}
	MATRICE_GLOBAL_FINL iterator_base& operator-=(difference_type _Offset) {
		return (*this += -(_Offset * _Mystep));
	}
	MATRICE_GLOBAL_FINL iterator_base operator-(difference_type _Offset) const {
		auto _Tmp = *this;
		return (_Tmp -= (_Offset));
	}
	MATRICE_GLOBAL_FINL iterator_base operator-(const iterator_base& _Right) const {
		return (_Myptr - _Right._Myptr);
	}
	MATRICE_GLOBAL_FINL reference operator[](difference_type _Offset) const {
		return (*(*this + _Offset * _Mystep));
	}
	MATRICE_GLOBAL_FINL bool operator==(const iterator_base& _Right) const {
		return (_Myptr == _Right._Myptr);
	}
	MATRICE_GLOBAL_FINL bool operator!=(const iterator_base& _Right) const {
		return (!(*this == _Right));
	}
	MATRICE_GLOBAL_FINL bool operator<(const iterator_base& _Right) const {
		return (_Myptr < _Right._Myptr);
	}
	MATRICE_GLOBAL_FINL bool operator>(const iterator_base& _Right) const {
		return (_Right < *this);
	}
	MATRICE_GLOBAL_FINL bool operator<=(const iterator_base& _Right) const {
		return (!(_Right < *this));
	}
	MATRICE_GLOBAL_FINL bool operator>=(const iterator_base& _Right) const {
		return (!(_Right > *this));
	}

	// \test for iterator end condition
	MATRICE_GLOBAL_FINL operator bool() { return (_Myptr != _Myend); }

	// \return pointer to current object
	MATRICE_GLOBAL_FINL operator pointer() { return (_Myptr); }

	// \forward range iteration methods for [this->_Myptr, this->_Myend)
	MATRICE_GLOBAL_FINL auto& begin() { return (*this); }
	MATRICE_GLOBAL_FINL auto end() { auto _Tmp = *this; return (_Tmp += _Mysize); }
	MATRICE_GLOBAL_FINL const auto& begin() const { return (*this); }
	MATRICE_GLOBAL_FINL const auto end() const { auto _Tmp = *this; return (_Tmp += _Mysize); }
	// \current iterator position
	MATRICE_GLOBAL_FINL auto pos() const { return std::distance(_Mybegin, _Myptr)/_Mystep; }

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
MATRICE_GLOBAL_FINL iterator_base<_Ty> operator+ (typename iterator_base<_Ty>::difference_type _Offset, iterator_base<_Ty> _Next) {
	return (_Next += _Offset);
}
/**********************************************************************
							Forward Range Iterator
	    Copyright (c) : Zhilong (Dgelom) Su, since 12/Jul/2018
 **********************************************************************/
template<typename _Ty, typename = std::enable_if_t<is_scalar_v<_Ty>>>
class _Matrix_forward_iterator : public iterator_base<_Ty>
{
	using _Mybase = iterator_base<_Ty>;
public:
	template<typename... _Args> 
	MATRICE_GLOBAL_FINL _Matrix_forward_iterator(_Args... _args)
		: _Mybase(_args...) {}

private:
	typename _Mybase::difference_type _Myidx, _Myidy;
};
/**********************************************************************
							Row-wise Forward Range Iterator
	    Copyright (c) : Zhilong (Dgelom) Su, since 12/Jul/2018
 **********************************************************************/
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
class _Matrix_rwise_iterator : public iterator_base<_Ty>
{
	/*<Note> 
		_Myptr  := input pointer to row
		_Mysize := number of rows
		_Mystep := number of cols
		_Return  := i-th element in the row that _Myptr point to
	  </Note>*/
	using base_type = iterator_base<_Ty>;
public:
	template<typename... _Args>
	MATRICE_GLOBAL_FINL _Matrix_rwise_iterator(_Args... _args) noexcept
		: base_type(_args...) {}

	// \element-wise accessor
	MATRICE_GLOBAL_FINL auto& operator[](typename base_type::difference_type _Offset) const {
		return (*(base_type::_Myptr + _Offset));
	}
	// \copy data from another iterator
	MATRICE_GLOBAL_FINL _Matrix_rwise_iterator& operator= (const base_type& _Iter) {
		std::copy(_Iter.begin(), _Iter.end(), base_type::_Myptr);
		return (*this);
	}
	// \copy data from another pointer
	MATRICE_GLOBAL_FINL _Matrix_rwise_iterator& operator= (const typename base_type::pointer _Ptr) {
		std::copy(_Ptr, _Ptr + base_type::_Mystep, base_type::_Myptr);
		return (*this);
	}
	// \copy data from initializer list
	MATRICE_GLOBAL_FINL _Matrix_rwise_iterator& operator= (const std::initializer_list<_Ty> _List) {
		std::copy(_List.begin(), _List.end(), base_type::_Myptr);
		return (*this);
	}
	// \iterator to tranverse the element in current row 
	MATRICE_GLOBAL_FINL auto begin() {
		return base_type(base_type::_Myptr, _Myrange);
	}
	MATRICE_GLOBAL_FINL const auto begin() const {
		return base_type(base_type::_Myptr, _Myrange);
	}
	MATRICE_GLOBAL_FINL auto end() {
		return base_type(base_type::_Myptr + _Myrange, _Myrange);
	}
	MATRICE_GLOBAL_FINL const auto end() const {
		return base_type(base_type::_Myptr + _Myrange, _Myrange);
	}
private:
	size_t _Myrange = base_type::_Mystep;
};

/**********************************************************************
						Column-wise Forward Range Iterator 
	    Copyright (c) : Zhilong (Dgelom) Su, since 12/Jul/2018
 **********************************************************************/
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
class _Matrix_cwise_iterator : public iterator_base<_Ty>
{
	/*<Note> 
		_Myptr  := input pointer to column
		_Mysize := number of cols
		_Mystep := 1
		_Return  := i-th element in the column that _Myptr point to
	  </Note>*/
	using base_type = iterator_base<_Ty>;
public:
	template<typename... _Args>
	MATRICE_GLOBAL_FINL _Matrix_cwise_iterator(_Args... _args) noexcept
		: base_type(_args...) {
		std::swap(_Myrange, base_type::_Mystep);
	}
	// \element-wise accessor
	MATRICE_GLOBAL_FINL auto& operator[](typename base_type::difference_type _Offset) const {
		return (*(base_type::_Myptr + _Offset*base_type::_Mysize));
	}
	// \copy data from another iterator
	MATRICE_GLOBAL_FINL _Matrix_cwise_iterator& operator= (const base_type& _Iter) {
		base_type _Tmp(base_type::_Myptr, _Myrange, base_type::_Mysize);
		std::copy(_Iter.begin(), _Iter.end(), _Tmp.begin());
		return (*this);
	}
	// \copy data from another pointer
	MATRICE_GLOBAL_FINL _Matrix_cwise_iterator& operator= (const typename base_type::pointer _Ptr) {
		base_type _Tmp(base_type::_Myptr, _Myrange, base_type::_Mysize);
		std::copy(_Ptr, _Ptr + base_type::_Mystep, _Tmp.begin());
		return (*this);
	}
	// \copy data from initializer list
	MATRICE_GLOBAL_FINL _Matrix_cwise_iterator& operator= (const std::initializer_list<_Ty> _List) {
		base_type _Tmp(base_type::_Myptr, _Myrange, base_type::_Mysize);
		std::copy(_List.begin(), _List.end(), _Tmp.begin());
		return (*this);
	}
	// \iterator to tranverse the element in current column
	MATRICE_GLOBAL_FINL auto begin() {
		return base_type(base_type::_Myptr, _Myrange, base_type::_Mysize);
	}
	MATRICE_GLOBAL_FINL const auto begin() const {
		return base_type(base_type::_Myptr, _Myrange, base_type::_Mysize);
	}
	MATRICE_GLOBAL_FINL auto end() {
		return base_type(_End(base_type::_Myptr, _Myrange, base_type::_Mysize), _Myrange, base_type::_Mysize);
	}
	MATRICE_GLOBAL_FINL const auto end() const {
		return base_type(_End(base_type::_Myptr, _Myrange, base_type::_Mysize), _Myrange, base_type::_Mysize);
	}
private:
	size_t _Myrange = 1;
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
