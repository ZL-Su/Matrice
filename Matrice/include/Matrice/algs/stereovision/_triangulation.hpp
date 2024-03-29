/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/
#pragma once
#include "core.hpp"
#include "solver.hpp"
#include "../geometry/transform.h"

MATRICE_ALGS_BEGIN
_DETAIL_BEGIN
template<typename _Ty>
class _Trig_base {
	static_assert(is_floating_point_v<_Ty>, 
		"_Ty in _Trig_base<_Ty> must be a floating point type.");
	using _Myt = _Trig_base;
public:
	using value_type = _Ty;
	using vector_type = Vec3_<value_type>;
	using matrix_2x8t = Matrix_<value_type, 2, 8>;
	using matrix_3x3t = Matrix_<value_type, 3, 3>;

	/**
	 * \brief Set relative geometry between two cameras.
	 * \param 'ext' a 6-vector holds the external parameters as {rx,ry,rz,tx,ty,tz}.
	 * \return the mutable reference of this object. 
	 */
	MATRICE_HOST_INL _Myt& set_relat_geo(initlist<value_type> ext) noexcept {
		decltype(auto) _It = ext.begin();
		const auto rx = *_It, ry = *(_It + 1), rz = *(_It + 2);
		const auto tx = *(_It + 3), ty = *(_It + 4), tz = *(_It + 5);

		rodrigues(vector_type{ rx, ry, rz }, m_rot[1]);
		m_trs[1].x = tx;
		m_trs[1].y = ty;
		m_trs[1].z = tz;
		
		return (*this);
	}
	/**
	 * \brief Set relative geometry between two cameras.
	 * \param 'r' a 3-vector holds the rotation vector {rx,ry,rz}.
	 * \param 't' a 3-vector holds the translation vector {tx,ty,tz}.
	 * \return the mutable reference of this object.
	 */
	MATRICE_HOST_INL _Myt& set_relat_geo(const vector_type& r, const vector_type& t) noexcept {
		rodrigues(r, m_rot[1]);
		m_trs[1].x = t.x;
		m_trs[1].y = t.y;
		m_trs[1].z = t.z;
		return (*this);
	}

	/**
	 * \brief Set internal geometry for each of the two. cams.
	 * \param 'kk' passes the internal parameters with a nested initializer list
	 * {{fx, fy, cx, cy}, {fx, fy, cx, cy}}.
	 * \return the mutable reference of this object.
	 */
	MATRICE_HOST_INL _Myt& set_inter_params(nested_initlist<value_type> kk) noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(kk.size() > 0 && kk.size() < 3, 
			"At least one set of internal parameters are required.");
#endif
		if (kk.size() == 1) {
			m_inpars.rview(0) = *kk.begin();
			m_inpars.rview(1) = m_inpars.rview(0);
		}
		else {
			m_inpars.rview(0) = *kk.begin();
			m_inpars.rview(1) = *(kk.begin() + 1);
		}
		return (*this);
	}

	/**
	 * \brief Remove image point distortion in in-place.
	 * \param 'i' the index of camera (view).
	 * \param 'u', 'v' the raw image point coordinates.
	 * \return the undistorted coordinates 'u', 'v'.
	 */
	MATRICE_HOST_INL void remove_distortion(size_t i, value_type& u, value_type& v) noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(i < m_inpars.rows(), "The index parameter 'i' overs the row range.");
#endif 
		const auto k = m_inpars[i];
		const auto fx = k[0], fy = k[1], cx = k[2], cy = k[3];
		const auto k1 = k[4], k2 = k[5], p1 = k[6], p2 = k[7];
		const auto x = (u - cx) / fx, y = (v - cy) / fy;
		const auto r2 = sq(x) + sq(y);
		const auto ud = 1 + k1 * r2 + k2 * sq(r2) + 2*(p1 * x + p2 * y);
		u = (ud * x + p2 * r2) * fx + cx;
		v = (ud * y + p1 * r2) * fy + cy;
	}

protected:
	matrix_2x8t m_inpars;
	matrix_3x3t m_rot[2];
	vector_type m_trs[2];
};

template<typename _Ty>
class _LSTrig : public _Trig_base<_Ty>
{
	using _Mybase = _Trig_base<_Ty>;
	using _Myt = _LSTrig;
public:
	using value_type = typename _Mybase::value_type;
	using vector_type = typename _Mybase::vector_type;
	using correspondence = pair_t<Vec2_<value_type>>;
	using correspondences = std::vector<correspondence>;
	_LSTrig() noexcept {
		_Mybase::m_rot[0].identity();
		_Mybase::m_trs[0] = 0;
		_Mybase::m_rot[1].identity();
		_Mybase::m_trs[1] = 0;
	}

	/**
	 * \brief set reference geometry between frames of the object and the ref. cam.
	 */
	MATRICE_HOST_INL _Myt& set_refer_geo(initlist<value_type> rt) noexcept {
		decltype(auto) _It = rt.begin();
		const auto rx = *_It, ry = *(_It + 1), rz = *(_It + 2);
		const auto tx = *(_It + 3), ty = *(_It + 4), tz = *(_It + 5);

		rodrigues(vector_type{ rx, ry, rz }, _Mybase::m_rot[0]);
		_Mybase::m_trs[0].x = tx;
		_Mybase::m_trs[0].y = ty;
		_Mybase::m_trs[0].z = tz;

		//compute the geometry bewtween the object and right cam frames.
		_Mybase::m_trs[1] = this->m_trs[1] + this->m_rot[1].mul(this->m_trs[0]);
		const auto r_temp = _Mybase::m_rot[1].mul(_Mybase::m_rot[0]).eval();
		_Mybase::m_rot[1] = r_temp;

		return (*this);
	}

	/**
	 * \brief Compute 3D coords with given stereo correspondence:
	 //tex:$\{\mathbf{p}\leftrightarrow \mathbf{q}\}$
	 * \param "p" and "q" the input point correspondence.
	 * \return the 3d point associated to the pair "p <--> q".
	 */
	MATRICE_HOST_INL auto compute(const Vec2_<value_type>& p, const Vec2_<value_type>& q) noexcept {
		Matrix_<value_type, compile_time_size<>::_4, compile_time_size<>::_3> A; 
		Vec4_<value_type> b;

		auto fx = _Mybase::m_inpars[0][0], fy = _Mybase::m_inpars[0][1];
		auto cx = _Mybase::m_inpars[0][2], cy = _Mybase::m_inpars[0][3];
		_Mybase::remove_distortion(0, p.x, p.y);
		auto x = p.x - cx, y = p.y - cy;
		auto R = _Mybase::m_rot[0];
		A.rview(0) = { x*R[2][0] - fx*R[0][0], 
			x*R[2][1] - fx*R[0][1], 
			x*R[2][2] - fx*R[0][2] };
		A.rview(1) = { y*R[2][0] - fy*R[1][0], 
			y*R[2][1] - fy*R[1][1], 
			y*R[2][2] - fy*R[1][2] };
		b(0) = fx * _Mybase::m_trs[0].x - x * _Mybase::m_trs[0].z;
		b(1) = fy * _Mybase::m_trs[0].y - y * _Mybase::m_trs[0].z;

		fx = _Mybase::m_inpars[1][0], fy = _Mybase::m_inpars[1][1];
		cx = _Mybase::m_inpars[1][2], cy = _Mybase::m_inpars[1][3];
		_Mybase::remove_distortion(1, q.x, q.y);
		x = q.x - cx, y = q.y - cy;
		R = _Mybase::m_rot[1];
		A.rview(2) = { x*R[2][0] - fx*R[0][0], 
			x*R[2][1] - fx*R[0][1], 
			x*R[2][2] - fx*R[0][2] };
		A.rview(3) = { y*R[2][0] - fy*R[1][0], 
			y*R[2][1] - fy*R[1][1], 
			y*R[2][2] - fy*R[1][2] };
		b(2) = fx * _Mybase::m_trs[1].x - x * _Mybase::m_trs[1].z;
		b(3) = fy * _Mybase::m_trs[1].y - y * _Mybase::m_trs[1].z;

		const auto solver = make_linear_solver<lls>(A);
		auto X = solver.solve(b);
		X(1) *= -1, X(2) *= -1;
		return (X);
	}
	MATRICE_HOST_INL auto compute(const pair_t<Vec2_<value_type>>& cpd) noexcept {
		return compute(cpd.first, cpd.second);
	}
};
_DETAIL_END
MATRICE_ALGS_END

DGE_MATRICE_BEGIN
enum trig_alg {
	LST = 0,
};
template<trig_alg _Alg, typename _Ty = default_type>
struct triangulation_engine {
	using type = conditional_t<is_same_constval_v<size_t(_Alg), size_t(trig_alg::LST)>, algs::detail::_LSTrig<_Ty>, void>;
};
template<trig_alg _Alg, typename _Ty = default_type>
using triangulation_t = typename triangulation_engine<_Alg, _Ty>::type;
DGE_MATRICE_END