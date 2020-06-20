/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once
#include <array>
#include "../util/_macros.h"
#include "../util/_std_wrapper.h"
#include "../util/_exception.h"
#include "../private/math/_primitive_funcs.hpp"
#include "_size_traits.h"

DGE_MATRICE_BEGIN
template<typename _Ity = size_t>
class shape_ : public std::array<_Ity, 4>
{
	using _Mybase = std::array<_Ity, 4>;
public:
	using index_type = typename _Mybase::value_type;

	shape_(initlist<index_type> shape) noexcept {
	}

	const index_type& h() const noexcept {
		return (*this)[0];
	}
	index_type& h() noexcept {
		return (*this)[0];
	}
	const index_type& w() const noexcept {
		return (*this)[1];
	}
	index_type& w() noexcept {
		return (*this)[1];
	}
	const index_type& d() const noexcept {
		return (*this)[2];
	}
	index_type& d() noexcept {
		return (*this)[2];
	}
};

template<size_t _Dim> struct shape_t {
	static_assert(_Dim==0||_Dim>3, 
		"_Dim in shape_t<_Dim> must be in the range [1, 3].");
};
template<> struct shape_t<1> {
	using value_type = size_t;

	constexpr auto size() noexcept { return h; }
	constexpr void reset() noexcept { h = 0; }
	constexpr auto rows() const noexcept { return h; }
	constexpr auto cols() const noexcept { return 1; }

	friend bool operator==(const shape_t<1>& _Left, const shape_t<1>& _Right) noexcept {
		return std::tie(_Left.h) == std::tie(_Right.h);
	}

	size_t h = 0; //height of the stored data
};
template<> struct shape_t<2> {
	using value_type = size_t;

	shape_t(size_t _h, size_t _w) noexcept
		:h(_h), w(_w) {
	}
	shape_t(const shape_t<2>& _other) noexcept
		:h(_other.h), w(_other.w) {
	}
	shape_t(shape_t<2>&& _other) noexcept
		:h(_other.h), w(_other.w) {
	}
	shape_t(initlist<size_t> _il) noexcept
		:h(*_il.begin()), w(*(_il.begin() + 1)) {
	}
	shape_t(shape_t<1>&& _other) noexcept
		:h(_other.h) {
	}

	constexpr auto size() noexcept { return h * w; }
	constexpr void reset() noexcept { h = w = 0; }
	constexpr auto rows() const noexcept { return h; }
	constexpr auto cols() const noexcept { return w; }
	/**
	 *\brief tile the 3d-shape to 2d-shape
	 */
	constexpr auto tile()const noexcept { return shape_t<1>{h*w}; }

	friend bool operator==(const shape_t<2>& _Left, const shape_t<2>& _Right) noexcept {
		return std::tie(_Left.h, _Left.w) == std::tie(_Right.h, _Right.w);
	}

	size_t h = 0; //height of the stored data
	size_t w = 0; //width of the stored data
};
template<> struct shape_t<3> {
	using value_type = size_t;
	using _Myt = shape_t<3>;
	shape_t(size_t _h, size_t _w, size_t _d=1) noexcept
		:h(_h), w(_w), d(_d) {
	}
	shape_t(const _Myt& _other) noexcept
		:h(_other.h), w(_other.w), d(_other.d) {
	}
	shape_t(_Myt&& _other) noexcept
		:h(_other.h), w(_other.w), d(_other.d){
	}
	shape_t(initlist<size_t> _il) noexcept
		:h(*_il.begin()), w(*(_il.begin()+1)), d(*(_il.begin() + 2)) {
	}
	shape_t(shape_t<2>&& _other) noexcept
		:h(_other.h), w(_other.w) {
	}
	shape_t(shape_t<1>&& _other) noexcept
		:h(_other.h) {
	}
	MATRICE_GLOBAL_FINL _Myt& operator=(const _Myt& _other) noexcept {
		h = _other.h, w = _other.w, d = _other.d;
		return (*this);
	}
	MATRICE_GLOBAL_FINL _Myt& operator=(_Myt&& _other) noexcept {
		h = _other.h, w = _other.w, d = _other.d;
		return (*this);
	}
	MATRICE_GLOBAL_FINL operator shape_t<2>() const noexcept {
		return shape_t<2>(h*d, w);
	}

	constexpr auto size() noexcept { return h * w * d; }
	constexpr void reset() noexcept { h = w = 0; d = 1; }
	constexpr auto rows() const noexcept { return h * d; }
	constexpr auto cols() const noexcept { return w; }
	/**
	 *\brief tile the 3d-shape to 2d-shape
	 */
	MATRICE_GLOBAL_FINL auto tile()const noexcept { 
		return shape_t<2>{h* d, w}; 
	}

	friend bool operator==(const _Myt& _Left, const _Myt& _Right) noexcept {
		return std::tie(_Left.h, _Left.w, _Left.d) == std::tie(_Right.h, _Right.w, _Right.d);
	}

	size_t h = 0; //height of the stored data
	size_t w = 0; //width of the stored data
	size_t d = 1; //depth of the stored data
};

_DETAIL_BEGIN
template<size_t _N>
MATRICE_GLOBAL_INL shape_t<_N> _Union(const shape_t<_N>& _1, const shape_t<_N>& _2) noexcept {
	if constexpr (_N == 1) 
		return shape_t<_N>{max(_1.h, _2.h)};
	if constexpr (_N == 2) 
		return shape_t<_N>{max(_1.h, _2.h), max(_1.w, _2.w)};
	if constexpr (_N == 3) 
		return shape_t<_N>{max(_1.h, _2.h), max(_1.w, _2.w), max(_1.d, _2.d)};
}
template<size_t _N1, size_t _N2>
MATRICE_GLOBAL_INL shape_t<max_integer_v<_N1, _N2>> _Union(const shape_t<_N1>& _1, const shape_t<_N2>& _2) noexcept {
	return _Union<max_integer_v<_N1, _N2>>(_1, _2);
}
_DETAIL_END
DGE_MATRICE_END
