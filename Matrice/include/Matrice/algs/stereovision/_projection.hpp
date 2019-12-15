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
#include "core.hpp"

MATRICE_ALGS_BEGIN
struct cs_alignment_tag { 
	struct left {}; struct middle {}; struct right{};
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
	}
	/**
	 *\brief initialize the rotation and translation
	 *\param [ext] external geometry parameters: rx, ry, rz, tx, ty, tz
	 */
	_Projection_base(const array_n<value_type, dim<<1>& ext) noexcept
		:m_tran(ext(3), ext(4), ext(5)) {
		point_type r{ ext(0), ext(1), ext(2) };
		const auto theta = sqrt(r.dot(r));
		const auto sval = sin(theta);
		auto cval = cos(theta);

		m_rotm.identity(); cval = 1 - cval; r = r / theta;

		m_rotm(0) += cval * r[0] * r[0];
		m_rotm(1) += cval * r[0] * r[1];
		m_rotm(2) += cval * r[0] * r[2];
		m_rotm(3) += cval * r[1] * r[0];
		m_rotm(4) += cval * r[1] * r[1];
		m_rotm(5) += cval * r[1] * r[2];
		m_rotm(6) += cval * r[2] * r[0];
		m_rotm(7) += cval * r[2] * r[1];
		m_rotm(8) += cval * r[2] * r[2];

		m_rotm(1) -= sval * r[2], m_rotm(2) += sval * r[1];
		m_rotm(3) += sval * r[2], m_rotm(5) -= sval * r[0];
		m_rotm(6) -= sval * r[1], m_rotm(7) += sval * r[0];

		m_irot = m_rotm.inv();
	}
	_Projection_base(const point_type& r, const point_type& t)noexcept
		: _Projection_base({r.x, r.y, r.z, t.x, t.y, t.z}) {
	}

	/**
	 *\brief rotation X with [x, y, z]^T = RX
	 *\param [_X] input 3d point
	 */
	MATRICE_HOST_INL point_type rotation(const point_type& _X)const noexcept{
		return (m_rotm.mul(_X) + m_tran);
	}
	/**
	 *\brief transform and normalize with [x, y, 1]^T = <RX + T>
	 *\param [_X] input 3d point
	 */
	MATRICE_HOST_INL point_type forward(const point_type& _X)const noexcept{
		point_type p = this->rotation(_X) + m_tran;
		return (p.normalize(p.z));
	}
	/**
	 *\brief back-projection
	 *\param [_x] input 2d image point [x, y, d]^T, where d is the depth.
	 */
	MATRICE_HOST_INL point_type backward(const point_type& _x)const noexcept {
		point_type& x = _x;
		x.x *= x.z, x.y *= x.z;
		x = x - m_tran;
		return m_irot.mul(x);
	}

protected:
	matrix_type m_rotm, m_irot;
	point_type m_tran;
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
	 *\brief compute the grad. of a re-projected point w.r.t. depth par.
	 *\param [pd] re-projected point, where the first two elements are coordinates of a normalized, distortion-rectified image point, the last one is its depth value.
	 */
	MATRICE_HOST_INL point_type grad(const point_type& pd)const noexcept {
		const auto X = _Mybase::rotation({ pd.x, pd.y, 1 });
		const auto s = 1 / sqr(X.z*pd.z + m_tran.z);
		const auto gx = (X.x*m_tran.z - X.z*m_tran.x)*s;
		const auto gy = (X.y*m_tran.z - X.z*m_tran.y)*s;
		return point_type{ gx, gy, 0 };
	}
	MATRICE_HOST_INL point_type grad(value_type x, value_type y, value_type depth)const noexcept {
		const auto X = _Mybase::rotation({ x, y, 1 });
		const auto s = 1 / sqr(X.z*depth + m_tran.z);
		const auto gx = (X.x*m_tran.z - X.z*m_tran.x)*s;
		const auto gy = (X.y*m_tran.z - X.z*m_tran.y)*s;
		return point_type{ gx, gy, 0 };
	}
private:
	using _Mybase::m_tran;
};
_DETAIL_END
template<typename _Ty>
using left_aligned_projection = detail::_Aligned_projection<_Ty, cs_alignment_tag::left>;
MATRICE_ALGS_END