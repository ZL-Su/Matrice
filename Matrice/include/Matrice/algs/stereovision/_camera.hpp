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

struct camera_tag {};

inline _DETAIL_BEGIN

template<typename _Ty, typename _Tag> class _Camera {
	static_assert(true, 
		"Invalid camera tag in dgelom::vision::detail::_Camera<_Ty, _Tag...>.");
};

/// <summary>
/// TRAIT, trait of a camera type
/// </summary>
template<typename _Ty, typename _Tag>
struct traits<_Camera<_Ty, _Tag>> {
	using value_type = _Ty;
	using category = _Tag;
	enum {
		// parameter size
		_Size = conditional_size_v<
		is_same_v<category, vision::persp_camera_tag>, 5,
		conditional_size_v<
		is_same_v<category, vision::ortho_camera_tag>, 4, 6>>
	};
};

/// <summary>
/// CLASS TEMPLATE, the base of camera types
/// </summary>
template<typename _Derived>
class _Camera<_Derived, camera_tag>  {
	using _Myt = _Camera;
	using _Mytraits = traits<_Derived>;
public:
	using value_type = typename _Mytraits::value_type;
	using category = typename _Mytraits::category;
	using init_list = std::initializer_list<value_type>;
	_Camera() = default;
	/**
	 * \brief CTOR, initialize camera with an internal calibration. 
	 */
	explicit _Camera(init_list _Calib) noexcept
		:_Mycalib(_Calib) {
	}
	explicit _Camera(init_list _Calib, init_list _Pose) noexcept
		:_Mycalib(_Calib), _Mypose(_Pose) {

	}

	/**
	 * \brief Eval forward projection. (Thread-safe method)
	 * \param 'X' 3D coodinates of an object point.
	 * \return auto [x, y, 1] = forward(X).
	 */
	MATRICE_HOST_INL auto forward(const Vec3_<value_type>& X) const noexcept {

	}

	MATRICE_HOST_INL auto backward(const Vec2_<value_type>& x) const noexcept {

	}

protected:
	static constexpr auto _Size = _Mytraits::_Size;

	// Camera calibration: $f_x, f_y, c_x, c_y$
	Vec_<value_type, _Size> _Mycalib;

	// Camera pose: $\r_x, r_y, r_z, t_x, t_y, t_z$
	Vec_<value_type, 6> _Mypose;

};

template<typename _Ty>
class _Camera<_Ty, persp_camera_tag> 
	: public _Camera<_Camera<_Ty, persp_camera_tag>, camera_tag>
{
	using _Myt = _Camera<_Ty, persp_camera_tag>;
	using _Mybase = _Camera<_Myt, camera_tag>;
	using _Mybase::_Size;
public:
	using typename _Mybase::value_type;

	_Camera() {};
	~_Camera() {};

private:

};

_DETAIL_END
MATRICE_ALG_END(vision)