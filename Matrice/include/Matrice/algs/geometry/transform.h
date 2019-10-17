/***********************************************************************
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

#include "core/matrix.h"
#include "core/vector.h"

DGE_MATRICE_BEGIN 
_DETAIL_BEGIN

template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
MATRICE_HOST_INL Vec3_<_Ty> _Rodrigues_impl(const Matrix_<_Ty, 3, 3>& _R) {
	using value_t = _Ty;
	constexpr const auto _Unit = value_t(1.0);
	constexpr const auto _Half = value_t(0.5);
	constexpr const auto _Zero = zero<value_t>::value;

	Vec3_<value_t> _Ret{_R[2][1]-_R[1][2], _R[0][2]-_R[2][0], _R[1][0]-_R[0][1]};

	auto s = _Ret.norm<2>() * _Half;
	auto c = (_R.trace() - _Unit)*_Half;
	c = c > _Unit ? _Unit : c < -_Unit ? -_Unit : c;
	auto a = std::acos(c);

	if (s < 1e-5) {
		if (c > _Zero) return (Vec3_<value_t>(_Zero));

		auto t = (_R[0][0] + _Unit)*_Half;
		_Ret.x = sqrt(max(t, _Zero));
		t = (_R[1][1] + _Unit)*_Half;
		_Ret.y = sqrt(max(t, _Zero))*(_R[0][1] < _Zero ? -_Unit : _Unit);
		t = (_R[2][2] + _Unit)*_Half;
		_Ret.z = sqrt(max(t, _Zero))*(_R[0][2] < _Zero ? -_Unit : _Unit);

		if (abs(_Ret.x)<abs(_Ret.y)&&abs(_Ret.x)<abs(_Ret.z)) {
			if ((_R[1][2] > _Zero) != (_Ret.y*_Ret.z > _Zero)) {
				_Ret.z = -_Ret.z;
			}
		}
		a /= _Ret.norm<2>();
		return (_Ret = a*_Ret);
	}
	
	return (_Ret = (_Half * a / s)*_Ret);
}

/**
 *\brief Compute rotation matrix between two 3d vectors
 *\param [_V1, _V2] the given two vectors
 */
template<typename _Ty>
MATRICE_HOST_INL Matrix_<_Ty, 3, 3> _Rot_from(const Vec3_<_Ty>& _V1, const Vec3_<_Ty>& _V2) {
	const auto _A = _V1.cross(_V2);
	const auto _R = Matrix_<_Ty, 3, 3>{0, _A[2],  _A[1], _A[2],     0, -_A[0],
		-_A[1], _A[0],     0
	};
	return decltype(_R)::diag(1) + _R + (_R*_R)*(1-_V1.dot(_V2)) / pow(_A.norm(), 2);
}

template<uint8_t, uint8_t> struct _Axis_type {};

/**
 *\brief TEMPLATE CLASS for axis-angle representation
 *\cite{https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation}
 */
template<typename, size_t> class _Axis_angle_rep {};

/**
 *\brief TEMPLATE CLASS for rigid transformation
 */
template<typename> class _GeoTransform_isometry {};

_DETAIL_END

template<typename _Ty, size_t _Cols, typename _Outty>
MATRICE_HOST_INL auto rodrigues(const Matrix_<_Ty, 3, _Cols>& _Left, _Outty _Right) {
	if constexpr (_Cols == 1) {

	}
	else {
		auto _Ret = detail::_Rodrigues_impl(_Left);
		if constexpr (std::is_pointer_v<_Outty>)
			dgelom::transform(_Ret.begin(), _Ret.end(), _Right);
		else _Right = _Ret;
	}
}

using axis_x_t = detail::_Axis_type<0, 3>;
using axis_y_t = detail::_Axis_type<1, 3>;
using axis_z_t = detail::_Axis_type<2, 3>;

template<typename _Ty, MATRICE_ENABLE_IF(is_floating_point_v<_Ty>)>
// *\brief TEMPLATE CLASS for rigid transformation
using isometry_t = detail::_GeoTransform_isometry<_Ty>;

template<typename _Ty, MATRICE_ENABLE_IF(is_floating_point_v<_Ty>)>
// *\brief TEMPLATE CLASS for axis-angle representation
using axisangle_t = detail::_Axis_angle_rep<_Ty, 3>;
DGE_MATRICE_END

#include "inline\_transform.inl"
