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
#include "../../../core"
DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty> class _Projection_base {
	using _Myt = _Projection_base<_Ty>;
public:
	static constexpr auto dim = 3;
	using value_type = _Ty;
	using point_type = Vec3_<value_type>;
	using matrix_type = Matrix_<value_type, dim, dim>;

	/**
	 *\brief initialize the rotation and translation
	 *\param [_Ext] external geometry parameters: rx, ry, rz, tx, ty, tz
	 */
	_Projection_base(const array_n<value_type, dim<<1>& _Ext)
	 :_MyT(_Ext(3), _Ext(4), _Ext(5)) {
		point_type r{ _Ext(0), _Ext(1), _Ext(2) };
		const auto theta = sqrt(r.dot(r));
		const auto sval = sin(theta);
		auto cval = cos(theta);

		_MyR.identity(); cval = 1 - cval; r = r / theta;

		_MyR(0) += cval * r[0] * r[0], _MyR(1) += cval * r[0] * r[1], _MyR(2) += cval * r[0] * r[2];
		_MyR(3) += cval * r[0] * r[1], _MyR(4) += cval * r[1] * r[1], _MyR(5) += cval * r[1] * r[2];
		_MyR(6) += cval * r[0] * r[2], _MyR(7) += cval * r[1] * r[2], _MyR(8) += cval * r[2] * r[2];

		_MyR(1) -= sval * r[2], _MyR(2) += sval * r[1];
		_MyR(3) += sval * r[2], _MyR(5) -= sval * r[0];
		_MyR(6) -= sval * r[1], _MyR(7) += sval * r[0];
	}

	/**
	 *\brief transform X with [x, y, z]^T = RX + T
	 *\param [_X] input 3d point
	 */
	MATRICE_HOST_INL point_type transform(const point_type& _X)noexcept{
		return (_MyR.mul(_X) + _MyT);
	}
	/**
	 *\brief transform and normalize with [x, y, 1]^T = <RX + T>
	 *\param [_X] input 3d point
	 */
	MATRICE_HOST_INL point_type forward(const point_type& _X)noexcept{
		point_type p = this->transform(_X);
		return (p.normalize(p.z));
	}

protected:
	matrix_type _MyR;
	point_type _MyT;
};
_DETAIL_END

DGE_MATRICE_END