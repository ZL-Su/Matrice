/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
***********************************************************************/
#pragma once
#include "../_matrix_base.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty>
struct index<_Ty, tag::_Tensor_tag> {
	static_assert(is_integral_v<_Ty>, "_Ty must be an integer type");

	MATRICE_GLOBAL_INL index(_Ty d, _Ty e, _Ty h, _Ty w) noexcept {
		data[0] = d, data[1] = e, data[2] = h, data[3] = w;
	}
	MATRICE_GLOBAL_INL _Ty value() noexcept {
		return 0;
	}
	_Ty data[4];
};
/**
 *\brief CLASS TEMPLATE dgelom::tensor prototype
 *\param <_Ty> data type
 *\param <_D> tensor depth which grows vertically
 *\param <_E> tensor extent in horizontal, the default is 1
 */
template<typename _Ty, size_t _D, size_t _E = 1>
class _Tensor_ : public types::Base_<_Tensor_<_Ty, _D, _E>, 
	tensor_traits<_Tensor_<_Ty, _D, _E>>>
{
	using _Myt = _Tensor_;
	using _Mytraits = tensor_traits<_Myt>;
	using _Mybase = types::Base_<_Myt, _Mytraits>;
public:
	enum { Size = 0, CompileTimeRows = 0, CompileTimeCols = 0 };
	using typename _Mybase::value_type;

	/**
	 *\brief constructor
	 *\param [_H, _W] rows and cols of each tensor cell
	 */
	_Tensor_(size_t _H, size_t _W) noexcept
		:_Mybase(_Mytraits::depth*_H, _Mytraits::extent*_W) {
		_Mybase::_Myshape = { _Mytraits::depth, _Mytraits::extent, _H, _W };
		_Mybase::_Flush_view_buf();
	}
	_Tensor_(size_t _H, size_t _W, value_type _Val) noexcept
		:_Tensor_(_H, _W) {
		_Mybase::operator=(_Val);
	}

	MATRICE_HOST_INL const value_type& operator()(size_t d, size_t e, size_t j, size_t i) const noexcept {
#ifdef _DEBUG || MATRICE_DEBUG
		DGELOM_CHECK(d < _D, "depth index over range.");
		DGELOM_CHECK(e < _E, "extent index over range.");
#endif
		const auto _Row = d * _Mybase::_Myshape.get(2) + j;
		const auto _Col = e * _Mybase::_Myshape.get(3) + i;
		return (_Mybase::operator[](_Row)[_Col]);
	}
	MATRICE_HOST_INL value_type& operator()(size_t d, size_t e, size_t j, size_t i) noexcept {
#ifdef _DEBUG || MATRICE_DEBUG
		DGELOM_CHECK(d < _D, "depth index over range.");
		DGELOM_CHECK(e < _E, "extent index over range.");
#endif
		const auto _Row = d * _Mybase::_Myshape.get(2) + j;
		const auto _Col = e * _Mybase::_Myshape.get(3) + i;
		return (_Mybase::operator[](_Row)[_Col]);
	}

	MATRICE_HOST_INL auto operator()(size_t d, size_t e) const noexcept {
#ifdef _DEBUG || MATRICE_DEBUG
		DGELOM_CHECK(d < _D, "depth index over range.");
		DGELOM_CHECK(e < _E, "extent index over range.");
#endif
		const auto w = _Mybase::_Myshape.get(3);
		const auto h = _Mybase::_Myshape.get(2);
		const auto r0 = d*h, r1 = (d+1) * h;
		const auto c0 = e*w, c1 = (e+1) * w;
		return _Mybase::block(c0, c1, r0, r1);
	}

private:
	size_t _Myw, _Myh;
};

template<typename _Ty, size_t _D, size_t _E>
struct tensor_traits<_Tensor_<_Ty, _D, _E>> {
	using type = _Ty;
	using category = tag::_Tensor_tag;
	static constexpr auto depth = _D;
	static constexpr auto extent = _E;
	static constexpr auto _M = 0, _N = 0;
	static constexpr bool Is_base = std::false_type::value;
};

_DETAIL_END
DGE_MATRICE_END