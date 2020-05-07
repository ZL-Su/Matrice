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
template<typename _Ty, size_t _Order> 
class _Spatial_transformer {};

template<typename _Ty>
class _Spatial_transformer<_Ty, 1> {
public:
	using value_type = _Ty;
	using param_type = Vec_<value_type, 4>;
	using jacob_type = Matrix_<value_type, 2, param_type::rows_at_compiletime>;
	struct warped_coords { value_type x, y; };

	/**
	 *\param [p] spatial deformation parameter
	 */
	_Spatial_transformer(const param_type& p) noexcept
		:_Mypar(p) {
	}

	/**
	 *\brief Transform local coords (x, y) to a new position.
	 */
	MATRICE_GLOBAL_INL auto operator()(value_type nx, value_type ny)const noexcept {
		return warped_coords {
			(1 + _Mypar[0])* nx + _Mypar[1] * ny, 
			_Mypar[2] * nx + (1 + _Mypar[3]) * ny
		};
	}

	/**
	 *\brief Update this transformer with:
	 //tex: $\mathbf{p}\leftarrow\mathbf{p}\circ\Delta\mathbf{p}^{-1}$
	 */
	MATRICE_GLOBAL_INL decltype(auto)update(const value_type* dp)noexcept {
		const auto _Det = safe_div(1, (1+dp[0])*(1+dp[3])-dp[1]*dp[2]);
		const auto idp11 = (1 + dp[3]) * _Det, idp12 = -dp[1] * _Det;
		const auto idp21 = -dp[2] * _Det, idp22 = (1 + dp[0]) * _Det;


		return (_Mypar);
	}

	/**
	 *\brief Get the Jacobian of this transformer w.r.t. the parameters
	 */
	MATRICE_GLOBAL_INL auto jacob(value_type x, value_type y) noexcept {
		return jacob_type{x, y, 0, 0, 0, 0, x, y};
	}

private:
	param_type _Mypar; //parameters are stored as $u_x, u_y, v_x, v_y$
};

_DETAIL_END
template<typename _Ty, size_t _Order>
using spatial_transformer_t = detail::_Spatial_transformer<_Ty, _Order>;
MATRICE_ALGS_END
