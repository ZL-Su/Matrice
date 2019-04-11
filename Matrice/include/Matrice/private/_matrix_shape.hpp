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

#include "../util/_macros.h"
#include "../util/_std_wrapper.h"
#include "../util/_exception.h"
#include "../private/math/_primitive_funcs.hpp"

DGE_MATRICE_BEGIN
/**
 * \2D shape type, auto [width, height] = shape(width, height)
 */
template<typename _Ity = size_t>
using shape_t = tuple<_Ity, _Ity>;
template<typename _Ity = size_t>
using shape3_t = tuple<_Ity, _Ity, _Ity>;
template<typename _Ity = size_t>
using shape4_t = tuple<_Ity, _Ity, _Ity, _Ity>;

/**
 *\brief CLASS TEMPLATE basic shape type
 *\param <_Ity> an integral template type
 */
template<typename _Ity = size_t,
	MATRICE_ENABLE_IF(is_integral_v<_Ity>)>
class basic_shape {
	using _Myt = basic_shape;
public:
	using value_type = _Ity;
	using view_shape = tuple<value_type,value_type>;
	using hist_shape = tuple<value_type,view_shape>;
	using full_shape = tuple<value_type,hist_shape>;

	MATRICE_GLOBAL_INL basic_shape(shape_t<value_type>&& _Shape) noexcept
		: _Data{ 1,{1,_Shape} } {}
	MATRICE_GLOBAL_INL basic_shape(shape3_t<value_type>&& _Shape) noexcept
		: _Data{ 1,{std::get<0>(_Shape), {std::get<1>(_Shape),std::get<2>(_Shape)}} } {}
	MATRICE_GLOBAL_INL basic_shape(shape4_t<value_type>&& _Shape) noexcept
		: _Data{ std::get<0>(_Shape), {std::get<1>(_Shape), {std::get<2>(_Shape),std::get<3>(_Shape)}} } {}
	template<typename _Jty>
	MATRICE_GLOBAL_INL basic_shape(std::initializer_list<_Jty> _Shape) noexcept {
		if (_Shape.size() == 2) {
			_Data = { 1,{1,{(value_type)*_Shape.begin(), (value_type)*(_Shape.begin() + 1)}} };
		}
		if (_Shape.size() == 3) {
			_Data = { 1,{(value_type)*_Shape.begin(),{(value_type)*(_Shape.begin() + 1), (value_type)*(_Shape.begin() + 2)}} };
		}
		if (_Shape.size() == 4) {
			_Data = { (value_type)*_Shape.begin(),{(value_type)*(_Shape.begin() + 1),{(value_type)*(_Shape.begin() + 2), (value_type)*(_Shape.begin() + 3)}} };
		}
	}
	MATRICE_GLOBAL_INL basic_shape(const _Myt& _Other) noexcept
		: _Data(_Other._Data) {}
	MATRICE_GLOBAL_INL basic_shape(_Myt&& _Other) noexcept
		: _Data(std::move(_Other._Data)) {}

	MATRICE_GLOBAL_INL auto& operator= (const _Myt& _Oth) {
		_Data = (_Oth._Data); return (*this);
	}
	MATRICE_GLOBAL_INL auto& operator= (_Myt&& _Oth) {
		_Data = std::move(_Oth._Data); return (*this);
	}
	MATRICE_GLOBAL_INL auto& operator= (const std::initializer_list<size_t> _Shape) {
		if (_Shape.size() == 2) {
			_Data = { 1,{1,{(value_type)*_Shape.begin(), (value_type)*(_Shape.begin() + 1)}} };
		}
		if (_Shape.size() == 3) {
			_Data = { 1,{(value_type)*_Shape.begin(),{(value_type)*(_Shape.begin() + 1), (value_type)*(_Shape.begin() + 2)}} };
		}
		if (_Shape.size() == 4) {
			_Data = { (value_type)*_Shape.begin(),{(value_type)*(_Shape.begin() + 1),{(value_type)*(_Shape.begin() + 2), (value_type)*(_Shape.begin() + 3)}} };
		}
		return (*this);
	}
	/**
	 *\brief Get unrolled shape data
	 */
	MATRICE_GLOBAL_INL auto operator()() const {
		return std::make_tuple(get(0), get(1), get(2), get(3));
	}
	/**
	 *\brief Get full rows and cols to a shape_t
	 */
	MATRICE_GLOBAL_INL operator shape_t<size_t>() const {
		return shape_t<size_t>(rows(), cols());
	}
	/**
	 *\brief Comparation operators
	 */
	MATRICE_GLOBAL_INL friend bool operator== (const _Myt& _lhs, const _Myt& _rhs) {
		return (_lhs._Data == _rhs._Data);
	}
	template<typename _Rhs>
	MATRICE_GLOBAL_INL friend bool operator== (const _Myt& _lhs, const _Rhs& _rhs) {
		return (_lhs._Data == _Myt(_rhs)._Data);
	}
	MATRICE_GLOBAL_INL friend bool operator!= (const _Myt& _lhs, const _Myt& _rhs) {
		return (_lhs._Data != _rhs._Data);
	}
	template<typename _Rhs>
	MATRICE_GLOBAL_INL friend bool operator!= (const _Myt& _lhs, const _Rhs& _rhs) {
		return (_lhs._Data != _Myt(_rhs)._Data);
	}
	/**
	 *\brief Get dim value at _Dim
	 */
	MATRICE_GLOBAL_INL constexpr auto get(uint8_t _Dim) const {
		if (_Dim == 0) return std::get<0>(_Data);
		if (_Dim == 1) return std::get<0>(std::get<1>(_Data));
		if (_Dim == 2) return std::get<0>(std::get<1>(std::get<1>(_Data)));
		if (_Dim == 3) return std::get<1>(std::get<1>(std::get<1>(_Data)));
		DGELOM_CHECK(_Dim<3, "_Dim over range of _Data.");
	}
	/**
	 *\brief Get full rows
	 */
	MATRICE_GLOBAL_INL constexpr auto rows() const {
		return get(0) * get(2);
	}
	/**
	 *\brief Get full cols
	 */
	MATRICE_GLOBAL_INL constexpr auto cols() const {
		return get(1) * get(3);
	}
	/**
	 *\brief return C*H*W
	 */
	MATRICE_GLOBAL_INL constexpr auto chw() const {
		return (get(1)*hw());
	}
	/**
	 *\brief return H*W
	 */
	MATRICE_GLOBAL_INL constexpr auto hw() const {
		return (get(2)*get(3));
	}
	/**
	 *\brief Parse the index for each dimension from a linear index
	 *\param [_Idx] input linear index
	 */
	MATRICE_GLOBAL_INL constexpr auto parse(size_t _Idx) const {
		auto n = _Idx / chw();  _Idx -= n * chw();
		auto c = _Idx / hw();   _Idx -= c * hw();
		auto h = _Idx / get(3); _Idx -= h * get(3);
		return shape4_t<>(n, c, h, _Idx);
	}

	friend MATRICE_GLOBAL_INL _Myt _Union(const _Myt& _1, const _Myt& _2) {
		return {max(_1.get(0), _2.get(0)), max(_1.get(1), _2.get(1)),
		max(_1.get(2), _2.get(2)), max(_1.get(3), _2.get(3)) };
	}

private:
	/**
	 *\brief formatted shape data: {N, {C, {H, W}}}
	 */
	full_shape _Data = { value_type(0),{value_type(0),{value_type(0),value_type(0)}} };
};

using basic_shape_t = basic_shape<>;

DGE_MATRICE_END

