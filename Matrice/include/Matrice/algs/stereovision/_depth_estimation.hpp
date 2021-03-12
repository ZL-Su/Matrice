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
#include "_spatial_transform.hpp"
#ifdef MATRICE_DEBUG
#include "../../io/io.hpp"
#endif

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
	_Depth_estimation_base(const initlist<value_t> rt)
		:m_projection{ rt } {
	}

	/**
	 *\brief getter/setter to the geometry of the stereovis. sys.
	 *\returns stereo projection geometry [m_projection]
	 */
	MATRICE_HOST_INL decltype(auto)geometry() const noexcept {
		return (m_projection);
	}
	MATRICE_HOST_INL decltype(auto)geometry() noexcept {
		return (m_projection);
	}

	/// <summary>
	/// Set to output verbose information or not.
	/// </summary>
	/// <returns>reference of verbose value [true or false]</returns>
	MATRICE_HOST_FINL decltype(auto)verbose()const noexcept {
		return m_verbose;
	}
	MATRICE_HOST_FINL decltype(auto)verbose() noexcept {
		return m_verbose;
	}

protected:
	MATRICE_HOST_INL void _Update_interpolator() noexcept {
		m_materp.reset(m_matching);
		m_referp.reset(m_reference);
	}

	/// <summary>
	/// Used to report some computation information.
	/// </summary>
	/// <param name="s">report string stream</param>
	MATRICE_HOST_INL void _Report(std::stringstream&& s) {
		/* 
		 s << "it=" << std::setfill('0') << std::setw(4) << iterations
		   << std::fixed << std::showpoint << std::setprecision(6)
		   << "    name1=" << val_1 << "    name2=" << val_2;
		 */
		MATRICE_USE_STD(cout);
		cout << s.str() << "\n";
	}

	projection_t m_projection;
	image_t m_reference, m_matching;
	interp_t m_referp, m_materp;

	bool m_verbose{ false };
};

template<typename _Ty, 
	class _Altag = cs_alignment_tag::left>
class _GCC_estimator 
	:public _Depth_estimation_base<_GCC_estimator<_Ty, _Altag>> {
	using _Myt = _GCC_estimator;
	using _Mybase = _Depth_estimation_base<_Myt>;
	using _Mystfer = _Spatial_transformer<_Ty, 1>;
	using _Mypart = Vec_<typename _Mybase::value_t, 4>;
	enum _Mytask {INIT_JACOB = 0, UPDATE_JACOB = 1};
public:
	// Define the order of spatial transformer
	static constexpr auto order = 1;
	// Define the size of the subsets
	static constexpr auto size = 41;
	using typename _Mybase::value_t;
	using typename _Mybase::point_t;
	using typename _Mybase::image_t;
	
	_GCC_estimator(const point_t& r, const point_t& t)
		:_Mybase{ r, t }, 
		m_jacob(sq<size_t>(size)), 
		m_resid(sq<size_t>(size)), 
		_Myle(size, size), 
		_Myri(size, size){
	}
	_GCC_estimator(const initlist<value_t> rt)
		:_Mybase{ rt }, 
		m_jacob(sq<size_t>(size)), 
		m_resid(sq<size_t>(size)), 
		_Myle(size, size), 
		_Myri(size, size){
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
	 *\param [x, y]: the coords of the image point in the left part of the stereo pair
	 */
	auto compute(value_t x, value_t y, value_t init_depth) {
#ifdef MATRICE_DEBUG
	    DGELOM_CHECK(!_Mybase::m_reference.empty, 
	    "_GCC_estimator<_Ty, _Altag> requires a valid stereo pair\n\
        be set with the method ::set_stereo_pair(l,r)\n\
        before computing the depth with ::compute(x, y, depth).");
#endif

		auto depth{ init_depth };
		_Mystfer T{ _Mypart::zeros() };
		
		//construct ref. subset surrounds (x, y)
		_Myle = _Mybase::m_reference.block(diff_t(x), diff_t(y), size/2);

		//compute the matching subset first time
		auto xp = _Mybase::m_projection.first(x, y, depth);
		_Myri = _Warp(T, xp.x, xp.y);
		array_n<value_t, ::dynamic> _Errmap(m_jacob.rows());
		_Errmap = _Myle - _Myri;

		//eval Jacobian and Hessian matrices
		m_jacob = _Diff<_Mytask::INIT_JACOB>(x, y, depth);
		auto hess = m_jacob.t().mul(m_jacob).eval();
		auto ihess = hess.spd().inv();

		//solve and update
		auto b = m_jacob.t().mul(_Errmap).eval();
		auto dp = ihess.mul(b).eval();
		dp = -1 * dp;
		auto residual = dp.norm<2>();
		depth += dp(0);
		T.update(dp.data() + 1);

		for (size_t k = 1; k < 50; ++k) {
			if (residual > 1.E-6) {
				xp = _Mybase::m_projection.update(depth);
#ifdef MATRICE_DEBUG
				print(xp);
#endif
				_Myri = _Warp(T, xp.x, xp.y);
				_Errmap = _Myle - _Myri;

				m_jacob = _Diff<_Mytask::UPDATE_JACOB>(x, y, depth);
				hess = m_jacob.t().mul(m_jacob);
				ihess = hess.spd().inv();
				
				b = m_jacob.t().mul(_Errmap);
				dp = ihess.mul(b);
				dp = -1 * dp;
				depth += dp(0);
				T.update(dp.data() + 1);
				residual = dp.norm<2>();
			}
		}

		return (depth);
	}

private:

	/**
	 *\brief Compute and update the Jacobian matrix.
	 *\param cx, cy: the x- and y-coordinates of the computation point.
	 *\param depth: the updated depth of (cx, cy).
	 */
	template<_Mytask _Task = _Mytask::INIT_JACOB>
	decltype(auto) _Diff(value_t cx, value_t cy, value_t depth) noexcept {
		const auto dxdp = _Mybase::m_projection.grad(cx, cy, depth);

		constexpr auto radius = size >> 1;
		for (diff_t dy = -radius; dy <= radius; ++dy) {
			const auto y = cy + dy;
			for (diff_t dx = -radius; dx <= radius; ++dx) {
				const auto x = cx + dx;
				const auto lidx = (dy + radius)*size + (dx + radius);
				auto J = m_jacob[lidx];
				if constexpr (_Task == _Mytask::UPDATE_JACOB) {
					const auto [dgdx, dgdy] = _Diff_g(x, y, depth);
					J[0] = -(dgdx * dxdp.x + dgdy * dxdp.y);
				}
				else if constexpr (_Task == _Mytask::INIT_JACOB) {
					const auto [dgdx, dgdy] = _Diff_g(x, y, depth);
					J[0] = -(dgdx * dxdp.x + dgdy * dxdp.y);

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
	value_t _Diff_d(value_t cx, value_t cy, value_t depth)const noexcept {
		const auto X = _Mybase::m_projection.backproj(cx, cy, depth);
		const auto p = _Mybase::m_projection.reproj(X);

		const auto [dgdx, dgdy] = _Mybase::m_materp.grad(p.x, p.y);
		const auto dxdp = _Mybase::m_projection.grad(cx, cy, depth);

		return -(dgdx*dxdp.x + dgdy * dxdp.y);
	}

	auto _Diff_g(value_t cx, value_t cy, value_t depth)const noexcept {
		const auto X = _Mybase::m_projection.backproj(cx, cy, depth);
		const auto p = _Mybase::m_projection.reproj(X);

		return _Mybase::m_materp.grad(p.x, p.y);
	}

	/**
	 *\brief Warp the subset:
	 //tex:$\Omega(\mathbf{x}'(d))$
	 */
	decltype(auto)_Warp(const _Mystfer& t, value_t xp, value_t yp)noexcept {
		const diff_t _Radius = _Myle.rows() >> 1;
		for (auto j = -_Radius; j <= _Radius; ++j) {
			auto ptr = _Myri[j + _Radius];
			for (auto i = -_Radius; i <= _Radius; ++i) {
				const auto [dx, dy] = t(i, j);
				ptr[i + _Radius] = _Mybase::m_materp(xp+dx, yp+dy);
			}
		}

		return (_Myri);
	}

private:
	// _Myle for the reference subset in the left image f 
	// _Myri is the moving subset in the right image g
	Matrix_<value_t, ::dynamic> _Myle, _Myri;
	/**
	 *\brief data arrangement:
	 //tex: $[\frac{\partial\varepsilon}{\partial d}, f_x\times dx, f_x\times dy, f_y\times dx, f_y\times dy]$
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
