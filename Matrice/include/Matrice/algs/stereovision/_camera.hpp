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
#include "_camera_pose.hpp"

MATRICE_ALG_BEGIN(vision)

/// <summary>
/// TAG, perspective camera
/// </summary>
struct persp_camera_tag {};
/// <summary>
/// TAG, pinhole camera
/// </summary>
struct pinho_camera_tag {};
/// <summary>
/// TAG, orthographic camera
/// </summary>
struct ortho_camera_tag {};
/// <summary>
/// TAG, refractive perspective camera
/// </summary>
struct refrap_camera_tag {};
/// <summary>
/// TAG, dummy camera
/// </summary>
struct camera_tag {};

MATRICE_HOST_INL
_DETAIL_BEGIN

// Forward declaration
template<typename _Ty, typename _Tag> 
class _Camera {
	static_assert(true, 
		"Unsupported camera tag '_Tag' in _Camera<_Ty, _Tag>.");
};

/// <summary>
/// TRAIT, trait of a camera type
/// </summary>
template<typename _Ty, typename _Tag>
struct traits<_Camera<_Ty, _Tag>> {
	using value_type = _Ty;
	using category = _Tag;
	enum {
		/// \brief *Number of parameters (fx, fy, cx, cy, [fs=0]).
		_Size = conditional_size_v<
		is_same_v<category, vision::persp_camera_tag>, 5,
		conditional_size_v<
		is_same_v<category, vision::ortho_camera_tag>, 4, 6>>
	};
};

template<class _Cam, typename _Ty = typename traits<_Cam>::value_type>
MATRICE_HOST_INL auto _Apply_distortion(const _Cam& cam, _Ty x, _Ty y) noexcept;
template<class _Cam, typename _Ty = typename traits<_Cam>::value_type>
MATRICE_HOST_INL auto _Remove_distortion(const _Cam& cam, _Ty u, _Ty v) noexcept;

/// <summary>
/// CLASS TEMPLATE, base of camera types.
/// </summary>
template<typename _Derived>
class _Camera<_Derived, camera_tag>  {
	using _Myt = _Camera;
	using _Mytraits = traits<_Derived>;
public:
	using value_type = typename _Mytraits::value_type;
	using category = typename _Mytraits::category;
	using init_list = MATRICE_STD(initializer_list)<value_type>;
	template<size_t N> using vector = auto_vector_t<value_type, N>;
	template<size_t N> using point = vector<N>;

	struct pose_type {
		MATRICE_HOST_INL pose_type(init_list pos) noexcept {
			r = pos.begin();
			t = pos.begin() + 3;
		}
		vector<3> r; // rotation vector (rx, ry, rz)
		vector<3> t; // translation vector (tx, ty, tz)

		MATRICE_GLOBAL_INL auto R() const noexcept {
			return rodrigues(r);
		}
		MATRICE_GLOBAL_INL auto matrix() const noexcept {
			return concat(this->R(), t);
		}
	};

	_Camera() = default;
	/**
	 * \brief CTOR, initialize camera with internal calibration. 
	 */
	explicit _Camera(init_list _Pars) noexcept
		:_Mypars(_Pars) {
	}
	/**
	 * \brief CTOR, initialize camera with internal calibration and pose.
	 */
	explicit _Camera(init_list _Pars, init_list _Pose) noexcept
		:_Camera(_Pars), _Mypose(_Pose) {
	}

	/**
	 * \brief Set or reset the camera pose.
	 * \param '_Pose' camera pose relative to the reference frame.
	 * \return Reference of the specified camera object.
	 */
	MATRICE_HOST_INL _Derived& pose(const pose_type& _Pose) noexcept {
		_Mypose = _Pose;
		return static_cast<_Derived>(*this);
	}
	MATRICE_HOST_INL decltype(auto) pose() const noexcept {
		return (_Mypose);
	}

	/**
	 * \brief Get projection matrix.
	 * \return A 3-by-4 projection matrix.
	 */
	MATRICE_HOST_INL auto pmatrix() noexcept {
		const auto _R = rodrigues(_Mypose.r);

	}

	/**
	 * \brief Eval forward projection to image domain (Thread-safe).
	 * \param 'X' 3D coodinates of an object point.
	 * \return auto [x, y, 1] = forward(X).
	 */
	MATRICE_HOST_INL auto forward(const point<3>& X) const noexcept {
		const auto _R = rodrigues(_Mypose.r);
		const auto _X = _R.mul(X) + _Mypose.t;
		return _X.eval<point<3>>();
	}
	/**
	 * \brief Eval backward projection to image domain (Thread-safe).
	 * \param 'p' observed pixel point.
	 * \return auto [x, y, 1] = backward(p).
	 */
	MATRICE_HOST_INL auto backward(const point<2>& p) const noexcept {
		const auto [x, y] = static_cast<const _Derived*>(this)->U(p);
		return point<3>{x, y, 1};
	}

protected:
	static constexpr auto _Size = _Mytraits::_Size;

	// Camera intrinsic parameters: $f_x, f_y, c_x, c_y [, f_s]$
	vector<_Size> _Mypars;

	// Camera pose with $r_x, r_y, r_z, t_x, t_y, t_z$
	pose_type _Mypose;

};

/// <summary>
/// \brief Specialization for perspective camera model
/// </summary>
/// <typeparam name="_Ty"></typeparam>
template<typename _Ty> 
requires is_floating_point_v<_Ty>
class _Camera<_Ty, persp_camera_tag> 
	: public _Camera<_Camera<_Ty, persp_camera_tag>, camera_tag> {
	using _Myt = _Camera<_Ty, persp_camera_tag>;
	using _Mybase = _Camera<_Myt, camera_tag>;
	using _Mybase::_Size;
public:
	using typename _Mybase::value_type;

	/**
	 * \brief GETER/SETTER, access to distortion model.
	 */
	decltype(auto) distortion_model() const noexcept {
		return _MyU;
	}
	decltype(auto) distortion_model() noexcept {
		return _MyU;
	}

	/**
	 * \brief Undistortion with distortion correction model 'U()'.
	 * \return [x, y] tuple of undistorted image coordinates in image domain.  
	 */
	MATRICE_HOST_INL auto U(value_type u, value_type v)const noexcept {
		const auto& _Myk = _Mybase::_Mypars;

		// backward transform in image space
		const auto fx = _Myk[0], fy = _Myk[1];
		const auto cx = _Myk[2], cy = _Myk[3];
		const auto x_d = (u - cx) / fx, y_d = (v - cy) / fy;

		// compute distortion correcting factor
		const auto r_2 = sqsum(x_d, y_d);
		const auto k1 = _MyU[0], k2 = _MyU[1];
		const auto p1 = _MyU[2], p2 = _MyU[3];
		const auto tmp = 1 + k1*r_2 + k2*sq(r_2) + 2*(p1*x_d + p2*y_d);

		// eval corrected image coordinates
		const auto x = tmp * x_d + p2 * r_2;
		const auto y = tmp * y_d + p1 * r_2;

		return tuple{ x, y };
	}
	/**
	 * \brief Distort an ideal image point [x, y] with distortion model 'D()'.
	 * \return [u, v] tuple of distorted pixel coordinates.
	 */
	MATRICE_HOST_INL auto D(value_type x, value_type y)const noexcept {
		const auto& _Myk = _Mybase::_Mypars;

		// compute distortion factor
		const auto r_2 = sqsum(x, y);
		const auto k1 = _MyD[0], k2 = _MyD[1];
		const auto p1 = _MyD[2], p2 = _MyD[3];
		const auto tmp = 1 + k1*r_2 + k2*sq(r_2) + 2*(p1*x + p2*y);

		// eval distorted image coordinates
		const auto x_d = tmp * x + p2 * r_2;
		const auto y_d = tmp * y + p1 * r_2;

		// forward transform in image space
		const auto fx = _Myk[0], fy = _Myk[1];
		const auto cx = _Myk[2], cy = _Myk[3];
		const auto u = x_d * fx + cx, v = y_d * fy + cy;

		return tuple{ u, v };
	}
private:
	// Distortion coefs: $k_1, k_2, p_1, p_2$
	auto_vector_t<value_type, 4> _MyU;
	auto_vector_t<value_type, 4> _MyD;
};

_DETAIL_END
/// <summary>
/// \brief ALIAS, Generic camera template
/// </summary>
/// <typeparam name="_Ty"></typeparam>
/// <typeparam name="_Tag"></typeparam>
template<typename _Ty, typename _Tag = persp_camera_tag>
using camera_t = detail::_Camera<_Ty, _Tag>;

/// <summary>
/// \brief ALIAS, Perspective camera template
/// </summary>
/// <typeparam name="_Ty">float or double, default is double</typeparam>
template<typename _Ty = default_type>
using perspective_camera_t = camera_t<_Ty, persp_camera_tag>;
/// <summary>
/// \brief ALIAS, Pinhole camera template
/// </summary>
template<typename _Ty = default_type>
using pinhole_camera_t = perspective_camera_t<_Ty>;

///<summary>
/// \brief FUNCTION, get camera category according to the tag of a given camera. 
///</summary>
template<typename _Ty, class _Tag>
MATRICE_HOST_FINL auto get_camera_category(camera_t<_Ty, _Tag>) noexcept {
	if constexpr (is_same_v<_Tag, persp_camera_tag>) {
		return "perspective camera";
	}
	if constexpr (is_same_v<_Tag, pinho_camera_tag>) {
		return "pinhole camera";
	}
	if constexpr (is_same_v<_Tag, ortho_camera_tag>) {
		return "orthographic camera";
	}
	if constexpr (is_same_v<_Tag, refrap_camera_tag>) {
		return "refractive pinhole camera";
	}
}

MATRICE_ALG_END(vision)