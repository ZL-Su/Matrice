/***********************************************************************
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

#include "core/matrix.h"
#include "core/vector.h"
#include "_geo_utils.hpp"

// clang-format off
MATRICE_GLOBAL_FINL constexpr auto operator""_mm(long double _Val) noexcept {
	return (_Val);
}
MATRICE_GLOBAL_FINL constexpr auto operator""_cm(long double _Val) noexcept {
	return (_Val);
}
MATRICE_GLOBAL_FINL constexpr auto operator""_m(long double _Val) noexcept {
	return (_Val);
}
MATRICE_GLOBAL_FINL constexpr auto operator""_s(long double _Val) noexcept {
	return (_Val);
}
MATRICE_GLOBAL_FINL constexpr auto operator""_degs(long double _Val) noexcept {
	return (_Val);
}
MATRICE_GLOBAL_FINL constexpr auto operator""_rads(long double _Val) noexcept {
	return (_Val);
}
// clang-format on

DGE_MATRICE_BEGIN
_DETAIL_BEGIN

template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
MATRICE_HOST_INL Vec3_<_Ty> _Rodrigues_impl(const Matrix_<_Ty, 3>& _R) {
	using value_t = _Ty;
	constexpr const auto _Unit = value_t(1.0);
	constexpr const auto _Half = value_t(0.5);
	constexpr const auto _Zero = value_t(0);

	Vec3_<value_t> _Ret{_R[2][1]-_R[1][2], _R[0][2]-_R[2][0], _R[1][0]-_R[0][1]};

	const auto s = _Ret.norm<2>() * _Half;
	auto c = (_R.trace() - _Unit)*_Half;
	c = c > _Unit ? _Unit : c < -_Unit ? -_Unit : c;
	auto a = MATRICE_STD(acos)(c);

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

/// <summary>
/// \brief Rodrigues transform from a vector to a matrix in 3-dimensional.
/// </summary>
/// <typeparam name="_Ty"> scalar </typeparam>
/// <param name="'_r'"> roration vector </param>
/// <returns>rotation matrix </returns>
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
MATRICE_HOST_INL Matrix_<_Ty, 3> _Rodrigues_impl(const Vec_<_Ty, 3>& w) noexcept {
	using value_t = _Ty;
	auto _Ret = Matrix_<value_t, 3>::diag(1);

	const auto theta = w.norm();
	if (theta != 0) {
		auto u = (w / theta).eval();
		const auto [s, c] = sin_cos(theta);
		_Ret(0) *= c, _Ret(4) *= c, _Ret(8) *= c;

		_Ret[0][0] += (1 - c) * u(0) * u(0);
		_Ret[0][1] += (1 - c) * u(0) * u(1) - s * u(2);
		_Ret[0][2] += (1 - c) * u(0) * u(2) + s * u(1);
		_Ret[1][0] += (1 - c) * u(1) * u(0) + s * u(2);
		_Ret[1][1] += (1 - c) * u(1) * u(1);
		_Ret[1][2] += (1 - c) * u(1) * u(2) - s * u(0);
		_Ret[2][0] += (1 - c) * u(2) * u(0) - s * u(1);
		_Ret[2][1] += (1 - c) * u(2) * u(1) + s * u(0);
		_Ret[2][2] += (1 - c) * u(2) * u(2);
	}
	return (_Ret);
}

/**
 * \brief Compute rotation matrix between two 3d vectors
 * \param '_V1, _V2' A pair of 3-vectors.
 */
template<typename _Ty> MATRICE_DEVICE_INL
Matrix_<_Ty, 3> _Rotation_between(const Vec3_<_Ty>& _V1, const Vec3_<_Ty>& _V2) {
	const auto _A = _V1.cross(_V2);
	const auto _R = Matrix_<_Ty, 3, 3>{0, _A[2],  _A[1], _A[2], 0, -_A[0],
		-_A[1], _A[0], 0};

	return decltype(_R)::diag(1) + _R + (_R*_R)*(1-_V1.dot(_V2)) / pow(_A.norm(), 2);
}

/// <summary>
/// \brief ENUM, defines geometric transform tags.
/// </summary>
enum class _Geotf_tag {
	ISO2D = 00, // 1d isometry or rigid-body transform
	ISO3D = 01, // 2d isometry or rigid-body transform
	SIM2D = 10, // similarity
	SIM3D = 11, // similarity
	AFF2D = 20, // affine
	AFF3D = 21, // affine
	PRO2D = 30, // projective
	PRO3D = 31, // projective
};

template<uint8_t, uint8_t> struct _Axis_type {};

/**
 * \brief TEMPLATE CLASS for axis-angle representation
 * \cite{https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation}
 */
template<typename, size_t> class _Axis_angle_rep {};

/**
 * \brief TEMPLATE CLASS for rigid transformation
 */
template<typename> class _Geotf_isometry {};

/// <summary>
/// \brief TEMPLATE CLASS, forward declaration for geometric transform.
/// </summary>
template<typename _Ty, _Geotf_tag... _Tag> class _Geo_transform{};

_DETAIL_END

/**
 * \brief FUNCTION TEMPLATE, convert degree angle to radian.
 * \param 'deg' input value with degree unit.
 * \return the angle value with the unit of radian.
 */
template<typename _Ty>
MATRICE_GLOBAL_FINL constexpr _Ty radian(_Ty deg) noexcept {
	return pi<_Ty> * deg / _Ty(180._degs);
}

/**
 * \brief FUNCTION TEMPLATE, convert degree angle to radian.
 * \param 'degs' input 3-vector with degree unit.
 * \return the angle vector with the unit of radian.
 */
template<typename _Ty>
MATRICE_GLOBAL_FINL constexpr Vec3_<_Ty> radian(Vec3_<_Ty> degs) noexcept {
	return { radian(degs[0]), radian(degs[1]), radian(degs[2]) };
}

/**
 * \brief FUNCTION TEMPLATE, conversion between the rotation vector and matrix.
 * \param "_In" an input 3d rotation vector or 3x3 rotation matrix.
 * \param "_Out" return type of this function.
 */
template<typename _Input, typename _Output>
MATRICE_HOST_INL auto rodrigues(const _Input& _In, _Output& _Out) noexcept {
	const auto _Ret = detail::_Rodrigues_impl(_In);
	if constexpr (is_pointer_v<_Output>)
		dgelom::transform(_Ret.begin(), _Ret.end(), _Out);
	else _Out = move(_Ret);
}

/**
 * \brief FUNCTION TEMPLATE, conversion between the rotation vector and matrix.
 * \param "_In" an input 3d rotation vector or 3x3 rotation matrix.
 * \return Rotation matrix or vector according to the input "_In".
 */
template<typename _Input>
MATRICE_HOST_INL auto rodrigues(const _Input& _In) noexcept {
	return detail::_Rodrigues_impl(_In);
}

// clang-format off
/// <summary>
/// \brief Define X axis in 3d space.
/// </summary>
using axis_x_t = detail::_Axis_type<0, 3>;
/// <summary>
/// \brief Define Y axis in 3d space.
/// </summary>
using axis_y_t = detail::_Axis_type<1, 3>;
/// <summary>
/// \brief Define Z axis in 3d space.
/// </summary>
using axis_z_t = detail::_Axis_type<2, 3>;

/// <summary>
/// \brief ALIAS TEMPLATE for axis-angle representation
/// </summary>
/// <typeparam name="_Ty"></typeparam>
template<typename _Ty, MATRICE_ENABLE_IF(is_floating_point_v<_Ty>)>
using axisangle_t = detail::_Axis_angle_rep<_Ty, 3>;

/// <summary>
/// \brief ALIAS of geometric transform tag.
/// </summary>
using geotf_tag = detail::_Geotf_tag;

/// <summary>
/// \brief ALIAS TEMPLATE for rigid transformation.
/// </summary>
/// <typeparam name="_Ty"></typeparam>
template<typename _Ty, MATRICE_ENABLE_IF(is_floating_point_v<_Ty>)>
using Isometry_t = detail::_Geotf_isometry<_Ty>;

/// <summary>
/// \brief ALIAS TEMPLATE for 2d/3d Euclidean transformation.
/// </summary>
/// <typeparam name="_Ty">Floating-point scalar type</typeparam>
/// <valueparam name="_Dim">Dimentionality of 2 or 3</valueparam>
template<typename _Ty, dim_t _Dim = 3D, 
	MATRICE_ENABLE_IF(is_floating_point_v<_Ty>)>
using isometry_t = detail::_Geo_transform<_Ty, _Dim == 2D ? geotf_tag::ISO2D : geotf_tag::ISO3D>;
// clang-format on
DGE_MATRICE_END
#include "inline\_transform.hpp"
