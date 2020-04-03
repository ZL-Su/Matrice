/*********************************************************************
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
***********************************************************************/
#pragma once
#include "core.hpp"

MATRICE_ALGS_BEGIN
_DETAIL_BEGIN
template<typename _Ty, size_t _Order, bool _Perfect> 
class _Spatial_transformer {};

template<typename _Ty>
class _Spatial_transformer<_Ty, 1, false> {
public:
	constexpr auto order = 0;
	using value_type = _Ty;
	using param_type = Vec_<value_type, 4>;
	using jacob_type = Matrix_<value_type, 2, param_type::rows_at_compiletime>;

	_Spatial_transformer(const param_type& param) noexcept
		:_Mypar(param) {
	}

	MATRICE_GLOBAL_INL auto warp(value_type x, value_type y) noexcept {
		return Vec2_<value_type>
		{
			(1 + _Mypar[0])* x + _Mypar[1] * y, 
			_Mypar[2] * x + (1 + _Mypar[3] * y)
		};
	}
	MATRICE_GLOBAL_INL auto jacob(value_type x, value_type y) noexcept {
		return jacob_type{x, y, 0, 0, 0, 0, x, y};
	}

private:
	param_type _Mypar; //$u_x, u_y, v_x, v_y$
};

_DETAIL_END
template<typename _Ty, size_t _Order>
using spatial_transformer_t = detail::_Spatial_transformer<_Ty, _Order>;
MATRICE_ALGS_END
