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
#include "core/vector.h"

MATRICE_ALGS_BEGIN
_DETAIL_BEGIN
struct _Signed_distance_function {
	template<typename _Ty>
	MATRICE_HOST_INL static _Ty sphere(const Vec3_<_Ty>& p) noexcept {
		return (sqrt(p.dot(p)) - one<_Ty>);
	}
	template<typename _Ty>
	MATRICE_HOST_INL static _Ty intersect(_Ty d1, _Ty d2) noexcept {
		return (max(d1, d2));
	}
	template<typename _Ty>
	MATRICE_HOST_INL static _Ty unions(_Ty d1, _Ty d2) noexcept {
		return (min(d1, d2));
	}
	template<typename _Ty>
	MATRICE_HOST_INL static _Ty diff(_Ty d1, _Ty d2) noexcept {
		return (max(d1, -d2));
	}
};
_DETAIL_END
/// <summary>
/// CLASS for signed distance function, including:
/// 1. 'sphere' -- 3d sphere distance;
/// 2. 'intersect' -- max of two values;
/// 3. 'unions' -- min of two values;
/// 4. 'diff' -- max of two values x and -y.
/// </summary>
using distance_function = detail::_Signed_distance_function;
MATRICE_ALGS_END