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
		:m_projection{ r,t } {
	}

	/**
	 *\brief getter/setter to the geometry of the stereovis. sys.
	 */
	MATRICE_HOST_INL decltype(auto) geometry() const noexcept {
		return (m_projection);
	}
	MATRICE_HOST_INL decltype(auto) geometry() noexcept {
		return (m_projection);
	}

protected:
	MATRICE_HOST_INL void _Update_interpolator() noexcept {
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
	using _Mypart = Vec_<typename _Mybase::value_t, 5>;
	enum _Mytask {INIT_JACOB = 0, UPDATE_JACOB = 1};
public:
	// Define the order of spatial transformer
	static constexpr auto order = 1;
	// Define the size of the subsets
	static constexpr auto size = 15;
	using typename _Mybase::value_t;
	using typename _Mybase::point_t;
	using typename _Mybase::image_t;
	
	_GCC_estimator(const point_t& r, const point_t& t)
		:_Mybase{ r, t }, 
		m_jacob(sqr<size_t>(size)), m_resid(sqr<size_t>(size)), 
		m_ref(size, size), m_rep(size, size){
	}

	/**
	 *\brief Set a stereo pair to be processed
	 *\param [left, right] left and right images of the stereo pair
	 */
	_Myt& set_stereo_pair(const image_t& left, const image_t& right) noexcept {
		_Mybase::m_reference = left;
		_Mybase::m_matching = right;
		_Mybase::_Update_interpolator();

		return (*this);
	}

	/**
	 *\brief Compute the depth for a given image point
	 *\param [x, y] the coords of the image point in the left part of the stereo pair
	 */
	auto compute(value_t x, value_t y, value_t init_depth) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(!_Mybase::m_reference.empty, "Set a stereo pair with method ::set_stereo_pair(l,r) before computing the depth in _GCC_estimator<>::compute(x, y, depth).");
#endif
		_Mypart pars = { init_depth, 0, 0, 0, 0 };

		//eval Jacobian and Hessiam matrices
		m_jacob = this->_Diff<_Mytask::INIT_JACOB>(x, y, pars[0]);
		auto hess = m_jacob.t().mul(m_jacob).eval();
		auto ihess = hess.spd().inv();
		
		//construct ref. subset surrounds (x, y)
		m_ref = _Mybase::m_reference.block(diff_t(x), diff_t(y), size >> 1);

		//compute the matching subset first time

	}

private:
	/**
	 *\brief Compute and update the Jacobian matrix.
	 */
	template<_Mytask _Task = _Mytask::INIT_JACOB>
	decltype(auto) _Diff(value_t rx, value_t ry, value_t depth) noexcept {
		constexpr auto radius = size >> 1;
		for (diff_t dy = -radius, r = 1; dy <= radius; ++dy, ++r) {
			const auto y = ry + dy;
			for (diff_t dx = -radius, c = 1; dx <= radius; ++dx, ++c) {
				const auto x = rx + dx;
				const auto lidx = (dy + radius)*size + (dx + radius);
				auto J = m_jacob[lidx];
				if constexpr (_Task == _Mytask::UPDATE_JACOB) {
					J[0] = _Diff_d(x, y, depth);
				}
				else if constexpr (_Task == _Mytask::INIT_JACOB) {
					J[0] = _Diff_d(x, y, depth);

					const auto[dfdx, dfdy] = _Mybase::m_referp.grad(x, y);
					J[1] = dfdx * dx, J[2] = dfdx * dy;
					J[3] = dfdy * dx, J[4] = dfdy * dy;
				}
				else { 
					static_assert("_Task must be _Mytask::INIT_JACOB or _Mytask::UPDATE_JACOB in method _Diff<>()"); 
				}
			}
		}
		return (m_jacob);
	}
	/**
	 *\brief Compute derivative of the residual w.r.t. the depth, e.g. 
	 //tex:$\dfrac{\partial\varepsilon}{\partial d}$
	 */
	value_t _Diff_d(value_t x, value_t y, value_t depth)const noexcept {
		const auto X = _Mybase::m_projection.backproj(x, y, depth);
		const auto xp = _Mybase::m_projection.reproj(X);
		const auto [dgdx, dgdy] = _Mybase::m_materp.grad(xp.x, xp.y);

		const auto dxdp = _Mybase::m_projection.grad(x, y, depth);

		return -(dgdx*dxdp.x + dgdy * dxdp.y);
	}

	auto _Residual() {

	}

private:
	// m_ref for the reference subset in f 
	// m_rep is the moving subset in g
	Matrix_<value_t, ::dynamic> m_ref, m_rep;
	/**
	 *\brief data arrangement:
	 //tex: $[\dfrac{\partial\varepsilon}{\partial d}, f_x\times dx, f_x\times dy, f_y\times dx, f_y\times dy]$
	 */
	Matrix_<value_t, ::dynamic, (order<<2) + 1> m_jacob;
	Matrix_<value_t, ::dynamic, 1> m_resid;
	point_t m_refpt, m_reppt;
	value_t m_depth = zero<value_t>;
};
template<typename _Ty, class _Altag>
struct _Estimator_traits<_GCC_estimator<_Ty, _Altag>> {
	using value_type = _Ty;
	using alignment_tag = _Altag;
	using projection_type = _Aligned_projection<value_type, alignment_tag>;
};

_DETAIL_END
template<typename _Ty>
using left_aligned_gcc_t = detail::_GCC_estimator<_Ty, cs_alignment_tag::left>;
MATRICE_ALGS_END
