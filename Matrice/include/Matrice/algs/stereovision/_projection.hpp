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
#include "../geometry/transform.h"

MATRICE_ALGS_BEGIN
struct cs_alignment_tag { 
	struct left {}; 
	struct middle {}; 
	struct right{};
};
_DETAIL_BEGIN
template<typename _Ty> class _Projection_base {
	using _Myt = _Projection_base<_Ty>;
public:
	static constexpr auto dim = 3;
	using value_type = _Ty;
	using point_type = Vec3_<value_type>;
	using matrix_type = Matrix_<value_type, dim>;

	template<typename _An> struct depth_paramterized_projection {
		point_type point;
		MATRICE_HOST_INL point_type grad() noexcept {
			return point_type{};
		}
	};

	_Projection_base() noexcept 
		:m_tran(0) {
		m_ints.rview(0) = zero<value_type>;
		m_ints.rview(1) = m_ints(0);
	}
	/**
	 *\brief initialize the rotation and translation
	 *\param [ext] external geometry parameters: rx, ry, rz, tx, ty, tz
	 */
	_Projection_base(const array_n<value_type, dim<<1>& ext) noexcept
		:m_tran{ ext(3), ext(4), ext(5) } {
		//set internal paramters to default
		m_ints.rview(0) = zero<value_type>;

		//cvt given external params to the rot. mat. and trans vec.
		rodrigues(point_type{ ext(0), ext(1), ext(2) }, m_rotm);

		m_irot = m_rotm.inv();
	}
	_Projection_base(const point_type& r, const point_type& t)noexcept
		: _Projection_base({r.x, r.y, r.z, t.x, t.y, t.z}) {
	}

	/**
	 *\brief getter/setter to internal parameters
	 */
	MATRICE_HOST_INL decltype(auto)internal_params()const noexcept {
		return (m_ints);
	}
	MATRICE_HOST_INL decltype(auto)internal_params() noexcept {
		return (m_ints);
	}

	/**
	 *\brief Rotate $\mathbf{X}$ with $[x, y, z]^T = \mathbf{RX}$
	 *\param [_X] input 3D point
	 */
	MATRICE_HOST_INL point_type rotate(const point_type& _X)const noexcept{
		return (m_rotm.mul(_X));
	}
	/**
	 *\brief Project a 3D point onto the image domain with: 
			 $[x, y, 1]^T = \mathbf{K}'<\mathbf{RX + T}>$
	 *\param [_X] input 3D point
	 */
	MATRICE_HOST_INL point_type forward(const point_type& _X)const noexcept{
		point_type p = this->rotate(_X) + m_tran;
		return (p.normalize(p.z));
	}

	/**
	 *\brief Inversely maps a given image point to its normalized counterpart.
	 *\param [_x] input 2D image point $[x, y, d]^T$, where $d$ is the depth.
	 */
	MATRICE_HOST_INL point_type backward(const point_type& _pd)const noexcept {
		return this->backward(_pd.x, _pd.y, _pd.z);
	}
	MATRICE_HOST_INL auto backward(value_type _x, value_type _y, value_type d)const noexcept {
		const auto ptr = m_ints[0];
		const auto fx = ptr[0], fy = ptr[1], cx = ptr[2], cy = ptr[3];
		const auto x = (_x - cx) / fx * d, y = (_y - cy) / fy * d;
		return point_type(x, y, d);
	}

protected:
	matrix_type m_rotm, m_irot;
	point_type m_tran;
	Matrix_<value_type, 2, 4> m_ints;
	value_type m_imx, m_imy;
};

template<typename _Ty, typename _An> class _Aligned_projection {
	static_assert("Unsupported coordinate system alignment type.");
};

template<typename _Ty>
class _Aligned_projection<_Ty, cs_alignment_tag::left>
	: public _Projection_base<_Ty> {
	using _Myt = _Aligned_projection;
	using _Mybase = _Projection_base<_Ty>;
public:
	using _Mybase::_Projection_base;
	using typename _Mybase::point_type;
	using typename _Mybase::value_type;

	/**
	 *\brief First time compute the reprojected counterpart of the reference image point with an initial depth.
	 *\param [x, y] ref. image coords to be lifted to the 3D space.
	 *\param [d] initially estimated depth value.
	 */
	MATRICE_HOST_INL auto first(value_type x, value_type y, value_type d)noexcept {
		_Mybase::m_imx = x, _Mybase::m_imy = y;
		return reproj(backproj(_Mybase::m_imx, _Mybase::m_imy, d));
	}

	/**
	 *\brief Update the reprojection with current depth.
	 *\param [depth] Current estimated depth value.
	 */
	MATRICE_HOST_INL auto update(value_type depth) const noexcept {
		return reproj(backproj(_Mybase::m_imx, _Mybase::m_imy, depth));
	}

	/**
	 *\brief compute the grad. of a re-projected point w.r.t. depth par.
	 *\param [pd] re-projected point, where the first two elements are coordinates of a normalized, distortion-rectified image point, the last one is its depth value.
	 */
	MATRICE_HOST_INL point_type grad(const point_type& pd)const noexcept {
		return this->grad(pd.x, pd.y, pd.z);
	}
	/**
	 *\brief Compute the derivaties:
	 //tex: $\dfrac{\partial \mathbf{x}'}{\partial d} = \mathbf{K}'\dfrac{\partial<\mathbf{RX}(d) + \mathbf{T}>}{\partial d}$
	 *\sa Eq.(3) in paper GCCA.
	 */
	MATRICE_HOST_INL point_type grad(value_type x, value_type y, value_type depth)const noexcept {
		const auto Rx = _Mybase::rotate({ x, y, 1 });

		const auto Tx = m_tran.x, Ty = m_tran.y, Tz = m_tran.z;
		const auto s = safe_div(1, sq(Rx.z * depth + Tz));
		const auto _Gxd = (Rx.x * Tz - Tx * Rx.z) * s;
		const auto _Gyd = (Rx.y * Tz - Ty * Rx.z) * s;

		const auto ptr = _Mybase::m_ints[1];
		return point_type{ ptr[0] * _Gxd, ptr[1] * _Gyd, 0 };
	}

	/**
	 *\brief back-project a given image point with its depth
	 *\param [pd] image coords with a depth in form of $[x, y, d]^T$
	 */
	template<typename... _Args>
	MATRICE_HOST_INL auto backproj(_Args&&... pd)const noexcept {
		return _Mybase::backward(pd...);
	}

	/**
	 *\brief reproject a given 3D point to the matching camera
	 *\param [_X] 3D coords
	 */
	MATRICE_HOST_INL auto reproj(const point_type& _X)const noexcept {
		auto q = _Mybase::forward(_X);
		const auto ptr = _Mybase::m_ints[1];
		const auto fx = ptr[0], fy = ptr[1], cx = ptr[2], cy = ptr[3];
		q.x = q.x * fx + cx, q.y = q.y * fy + cy;
		return (q);
	}
private:
	using _Mybase::m_tran;
};
_DETAIL_END
template<typename _Ty>
using left_aligned_projection = detail::_Aligned_projection<_Ty, cs_alignment_tag::left>;
MATRICE_ALGS_END