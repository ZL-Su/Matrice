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

	using param_v6 = Vec_<value_type, 6>;
	using jacob_v6 = Matrix_<value_type, 2, 5>;

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
	 //tex: $\mathcal{T}(\eta;\mathbf{p})\leftarrow\mathcal{T}(\eta;\mathbf{p})\circ\mathcal{T}^{-1}(\eta;\Delta\mathbf{p})$
	 */
	MATRICE_GLOBAL_INL decltype(auto)update(const value_type* dp)noexcept {
		const auto idp0 = 1 + dp[3], idp1 = -dp[1];
		const auto idp2 = -dp[2], idp3 = 1 + dp[0];
		const auto _Det = safe_div(1, idp0*idp3-idp1*idp2);
		
		auto p = _Mypar.data();

		const auto p0p1 = 1 + p[0];
		const auto p0 = _Det*(p0p1 * idp0 + p[1] * idp2) - 1;
		const auto p1 = _Det*(p0p1 * idp1 + p[1] * idp3);

		const auto p3p1 = 1 + p[3];
		const auto p2 = _Det*(p[2] * idp0 + p3p1 * idp2);
		const auto p3 = _Det*(p[2] * idp1 + p3p1 * idp3) - 1;

		p[0] = p0, p[1] = p1, p[2] = p2, p[3] = p3;
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
