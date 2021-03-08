/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
*********************************************************************/
#pragma once

#include <core.hpp>

MATRICE_ALG_BEGIN(vision)

/// <summary>
/// TAG, perspective camera
/// </summary>
struct persp_camera_tag {};
/// <summary>
/// TAG, orthographic camera
/// </summary>
struct ortho_camera_tag {};
/// <summary>
/// TAG, refractive perspective camera
/// </summary>
struct refrap_camera_tag {};

inline _DETAIL_BEGIN
template<typename _Ty, typename _Tag>
requires is_floating_point_v<_Ty> class _Camera {
	using _Myt = _Camera;
public:
	using value_type = _Ty;
	using category = _Tag;
	
	/**
	 * \brief Eval forward projection.
	 * \param 'X' 3D coodinates of an object point.
	 * \return auto [x, y, 1] = forward(X).
	 */
	MATRICE_HOST_INL auto forward(const Vec3_<value_type>& X) const noexcept {

	}

protected:
	// Camera calibration: $f_x, f_y, c_x, c_y$
	Vec_<value_type, 4> _Mycalib;

	// Camera pose: $\r_x, r_y, r_z, t_x, t_y, t_z$
	Vec_<value_type, 6> _Mypose;

};
_DETAIL_END
MATRICE_ALG_END(vision)