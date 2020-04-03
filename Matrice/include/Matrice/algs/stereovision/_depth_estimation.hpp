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
#include "algs/correlation/_optim.h"
#include "algs/interpolation/_cubic_conv_interp.hpp"
#include "_projection.hpp"

MATRICE_ALGS_BEGIN
_DETAIL_BEGIN
template<typename _Est> struct _Estimator_traits {
	using value_type = default_type;
	using alignment_tag = cs_alignment_tag::left;
	using projection_type = _Aligned_projection<value_type, alignment_tag>;
};

template<typename _Derived>
class _Depth_estimation_base {
	using _Myt = _Depth_estimation_base;
	using _Mytraits = _Estimator_traits<_Derived>;
public:
	using value_t = typename _Mytraits::value_type;
	using projection_t = typename _Mytraits::projection_type;
	using point_t = typename projection_t::point_type;
	using interp_t = interpolation<value_t, bilerp_tag>;
	using image_t = Matrix<value_t>;

	_Depth_estimation_base(const point_t& r, const point_t& t)
		:m_projection(r,t){
	}

protected:
	void _Update_interpolator() noexcept {
		m_materp.set(m_matching);
		m_referp.set(m_reference);
	}
	projection_t m_projection;
	image_t m_reference, m_matching;
	interp_t m_referp, m_materp;
};

template<typename _Ty, 
	class _Altag = cs_alignment_tag::left>
class _GCC_estimator 
	:public _Depth_estimation_base<_GCC_estimator<_Ty, _Altag>> {
	using _Myt = _GCC_estimator;
	using _Mybase = _Depth_estimation_base<_Myt>;
	enum _Mytask {INIT_JACOB = 0, UPDATE_JACOB = 1};
public:
	static constexpr auto order = 1; // shape function order
	static constexpr auto size = 15; // subset size
	using typename _Mybase::value_t;
	using typename _Mybase::point_t;
	using typename _Mybase::image_t;
	
	_GCC_estimator()
		:_Mybase(point_t(), point_t()), m_jacob(sqr<size_t>(size)) {
	}
	_GCC_estimator(const point_t& r, const point_t& t)
		:_Mybase(r, t), m_jacob(sqr<size_t>(size)) {
	}

	_Myt& set_stereo_pair(const image_t& left, const image_t& right) noexcept {
		_Mybase::m_reference = left;
		_Mybase::m_matching = right;
		_Mybase::_Update_interpolator();

		return (*this);
	}

private:
	/**
	 *\brief computes and updates Jacobian matrix.
	 */
	template<_Mytask _Task = _Mytask::INIT_JACOB>
	void _Diff() noexcept {
		constexpr auto radius = size >> 1;
		for (diff_t dy = -radius; dy <= radius; ++dy) {
			const auto y = m_refpt.y + dy;
			for (diff_t dx = -radius; dx <= radius; ++dx) {
				const auto x = m_refpt.x + dx;
				const auto lidx = (dy + radius)*size + (dx + radius);
				if constexpr (_Task = _Mytask::UPDATE_JACOB) {
					m_jacob[lidx][6] = _Diff_d(x, y, m_depth);
				}
				else if constexpr (_Task = _Mytask::INIT_JACOB) {
					const auto[dfdx, dfdy] = _Mybase::m_referp.grad(x, y);
					auto J = m_jacob[lidx];
					J[0] = dfdx, J[1] = dfdx * dx, J[2] = dfdx * dy;
					J[3] = dfdy, J[4] = dfdy * dx, J[5] = dfdy * dy;
					J[6] = _Diff_d(x, y, m_depth);
				}
				else { 
					static_assert("_Task must be _Mytask::INIT_JACOB or _Mytask::UPDATE_JACOB in method _Diff<>()"); 
				}
			}
		}
	}
	/**
	 *\brief computes derivative of the residual w.r.t. the depth.
	 */
	value_t _Diff_d(value_t x, value_t y, value_t depth)const noexcept {
		const auto X = _Mybase::m_projection.backward({ x, y, depth });
		const auto xp = _Mybase::m_projection.forward(X);
		const auto dxdp = _Mybase::m_projection.grad(x, y, depth);
		const auto [dgdx, dgdy] = _Mybase::m_materp.grad(xp.x, xp.y);

		return -(dgdx*dxdp.x + dgdy * dxdp.y);
	}

private:

	/**
	 *\data arrangement:
	 * [ f_x f_x*dx f_x*dy f_y f_y*dx f_y*dy e_d ]
	 */
	Matrix_<value_t, ::dynamic, order + 1> m_jacob;
	point_t m_refpt, m_reppt;
	value_t m_depth = zero<value_t>;
};
template<typename _Ty, class _Altag>
struct _Estimator_traits<_Stereo_icgn_estimator<_Ty, _Altag>> {
	using value_type = _Ty;
	using alignment_tag = _Altag;
	using projection_type = _Aligned_projection<value_type, alignment_tag>;
};

_DETAIL_END
MATRICE_ALGS_END
