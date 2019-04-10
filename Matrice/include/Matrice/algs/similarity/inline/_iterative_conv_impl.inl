/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
#include "../_iterative_conv_impl.h"
#include "../../../private/_range.h"
#include "../../../private/nonfree/_lnalge.h"

MATRICE_ALGS_BEGIN namespace detail {
namespace conv_internal {
/**
 *\compute sequence: _Vec = {_Val*x^0, _Val*x^1, ..., _Val*x^n}
 */
template<typename _Vecty> MATRICE_HOST_FINL
void _Delta_pow_n(_Vecty& _Vec, typename _Vecty::value_t _Val) {
	for (auto _It = _Vec.begin() + 1; _It != _Vec.end(); ++_It) {
		*_It = _Val * *(_It - 1);
	}
}

/**
 * \determines region boundary from given center and radius
 * \return: [left, right, up, down] 
 */
template<typename _Pty, typename _Ty = typename _Pty::value_t>
MATRICE_HOST_INL auto _Square_range(const _Pty& _Center, const _Ty& _Radius) {
	int _L = floor<int>(_Center.x) - _Radius;
	int _R = floor<int>(_Center.x) + _Radius + 1;
	int _U = floor<int>(_Center.y) - _Radius;
	int _D = floor<int>(_Center.y) + _Radius + 1;

	return std::make_tuple(_L, _R, _U, _D);
}

/**
 *\interpolators
 */
template<size_t _N1, size_t _N2, typename... _Ts> struct _Op {
	template<typename... _Args> auto operator()(const _Args&... _Args) {}
};
template<typename _Maty, typename _Vecty>
struct _Op<bcspline, 1, _Maty, _Vecty> {
	MATRICE_HOST_FINL _Op(const _Maty& _M, const _Vecty& _V1, const _Vecty& _V2) : m_coef(_M), m_xbuf(_V1), m_ybuf(_V2) {}

	// \calc. interpolation value at (x, y) in current image
	MATRICE_HOST_FINL auto operator()(int x, int y) const {
		const auto _Row_0 = m_coef[y    ], _Row_1 = m_coef[y + 1];
		const auto _Row_2 = m_coef[y + 2], _Row_3 = m_coef[y + 3];
		const auto _Row_4 = m_coef[y + 4], _Row_5 = m_coef[y + 5];
		
		return (_Row_0[x] * m_xbuf[0] + _Row_1[x] * m_xbuf[1] + _Row_2[x]);
	}

	const _Maty& m_coef;
	const _Vecty& m_xbuf, m_ybuf;
};

/**
 *\estimate Jacobian: df/d\mathbf{p} = $\partial f / \partial p_i$
 */
template<> struct _Op<1, 6> {
	template<typename _Mulmty, typename _Poty, typename _Rgty>
	static auto J(const _Mulmty& _Ref, const _Poty& _Pos, const _Rgty& _Rx, const _Rgty& _Ry) {
		using value_type = typename _Mulmty::value_type;

		//auto [_Dfdx, _Dfdy] = _Ref.view_n<2,1>(_Rx, _Ry);
		const auto &_Dfdx = _Ref[1], &_Dfdy = _Ref[2];
		const auto _Size = max(_Rx.size(), _Ry.size());

		tensor<value_type, 1, 6> _Ret(_Size, _Size);
		for (auto y : _Ry) { auto _Dy = y - _Pos.y;
			const auto _y = y - _Ry.begin();
			for (auto x : _Rx) { auto _Dx = x - _Pos.x;
				const auto _x = x - _Rx.begin();
				const auto _Fx = _Dfdx[y][x], _Fy = _Dfdy[y][x];
				_Ret(_y, _x) = { _Fx, _Fx*_Dx, _Fx*_Dy, _Fy, _Fy*_Dx, _Fy*_Dy };
			}
		}
		return std::forward<decltype(_Ret)>(_Ret);
	};
};

template<std::size_t _ORDER> struct _Jacobian {
	static_assert(_ORDER > 2, "Oops, unsupported warp order: _ORDER.");
};

template<> struct _Jacobian<1> {
	template<typename _Mty, typename _Pty, typename _Rty>
	MATRICE_HOST_INL static auto eval(const _Mty& _Ref, const _Pty& _Pos, const _Rty& _Rx, const _Rty& _Ry) {
		using value_type = typename _Mty::value_type;
		const auto _Size = max(_Rx.size(), _Ry.size());
		auto[_dfdx, _dfdy] = _Ref.view_n<2, 1>(_Rx, _Ry);

		tensor<value_type, 1, 6> _Ret(_Size, _Size);
		for (auto y : _Ry) {
			const auto _Dy = y - _Pos.y;
			const auto _y = y - _Ry.begin();
			auto _Fx = _dfdx[y], _Fy = _dfdy[y];
			for (auto x : _Rx) {
				const auto _Dx = x - _Pos.x;
				const auto _x = x - _Rx.begin();
				_Ret(_y, _x) = { _Fx[x], _Fx[x]*_Dx, _Fx[x]*_Dy, 
					_Fy, _Fy[x]*_Dx, _Fy[x]*_Dy };
			}
		}
		return std::forward<decltype(_Ret)>(_Ret);
	};
};

template<std::size_t _ORDER> struct _Warp_update {};

template<> struct _Warp_update<compile_time_size<>::val_1> {
	template<typename _Vecty> 
	static MATRICE_HOST_INL auto& fwd(_Vecty& x, const _Vecty& y) {}

	/**
	 * \inverse compositional warp parameter update.
	 * \ x := p(u, dudx, dudy, v, dvdx, dvdy)
	 * \ y := update of parameter vector p, overwritten by the updated p
	 */
	template<typename _Vecty>
	static MATRICE_HOST_INL auto& inv(const _Vecty& x, _Vecty& y) {
		auto _Inv_det = 1/((1 + y[1])*(1 + y[5]) - y[2] * y[4]);
		auto _Val_0 = y[0] * (1 + y[5]) - y[3] * y[2];
		auto _Val_1 = y[3] * (1 + y[1]) - y[0] * y[4];
		auto _Val_2 = y[5] + 1, _Val_3 = y[4];
		auto _Val_4 = y[3] + 1, _Val_5 = y[2];

		y[0] = x[0] - ((1 + x[1])*_Val_0 - x[2] * _Val_1)*_Inv_det;
		y[1] = ((1 + x[1])*_Val_2 - x[2] * _Val_3)*_Inv_det - 1;
		y[2] = (x[2] * _Val_4 + (1 + x[1])*_Val_5)*_Inv_det;
		y[3] = x[3] - ((1 + x[5])*_Val_1 - x[4] * _Val_0)*_Inv_det;
		y[4] = (x[4] * _Val_2 - (1 + x[5])*_Val_3)*_Inv_det;
		y[5] = ((1 + x[5])*_Val_4 - x[4] * _Val_5)*_Inv_det - 1;

		return (y);
	}
};
}

/**
 * \update deformed subset for each interation, return 2/\sqrt{\sum{(g-\bar{g})^2}}
 */
template<typename _Derived> MATRICE_HOST_FINL 
auto _Iterative_conv_base<_Derived>::_Update_subset(const param_type& _P) {
#define _WITHIN_RANGE_OF_REFIMG \
	if (_Iy - _Bl >= 0 && _Ix - _Bl >= 0 && \
		 _Iy + _Bu < m_reference[0].rows() && _Ix + _Bu < m_reference[0].cols())

	constexpr auto _Bl = _Conv_border_size<interp_category>::lower;
	constexpr auto _Bu = _Conv_border_size<interp_category>::upper;
	const auto _Radius = static_cast<int>(m_options());
	const auto [_L, _R, _U, _D] = conv_internal::_Square_range(m_pos, _Radius);

	stack_vector _Buf_x{ 1. }, _Buf_y{ 1. };

	auto _Mean = zero<value_type>;
	for (auto y = _U, _Off_y = 0; y < _D; ++y, ++_Off_y) {
		auto _Dy = static_cast<value_type>(y - m_pos.y);

		auto _X = _L + _P[0] + _P[2] * _Dy;
		auto _Y =  y + _P[3] + _P[5] * _Dy;
		for (auto x = _L, _Off_x = 0; x < _R; ++x, ++_Off_x) {
			auto _Dx = static_cast<value_type>(x - m_pos.x);
			_X += 1 + _P[1] * _Dx; _Y += _P[4] * _Dx;

			auto _Ix = floor<int>(_X), _Iy = floor<int>(_Y);
			_WITHIN_RANGE_OF_REFIMG {
				//auto _Delta_x = _X - static_cast<value_type>(_Ix);
				//auto _Delta_y = _Y - static_cast<value_type>(_Iy);
				//conv_internal::_Delta_pow_n(_Buf_x, _Delta_x);
				//conv_internal::_Delta_pow_n(_Buf_y, _Delta_y);

				_Mean += _Mycur(_Off_y, _Off_x) = (*_Myitp)(_Ix-_Bl, _Iy);
			}
			else _Mycur(_Off_y, _Off_x) = zero<value_type>;
		}
	}
	_Mean /= static_cast<value_type>(multiply(_R - _L, _D - _U));

	_Mycur = _Mycur - _Mean;
	auto _SSD = sqrt((_Mycur*_Mycur).sum());
	if (_SSD < m_options)
		throw std::runtime_error("Bad value of variable _SSD in _Invcomp_conv_base<...>::_Update_subset(_P).");
	_Mycur = _Mycur / _SSD;

	return (2./_SSD);

#undef _WITHIN_RANGE_OF_REFIMG
}

/**
 * \For each point refinement, this method aims to compute the Jacobian and the Hessian before stepping into the iterative solver _Impl().
 */
template<typename _Ty, typename _Tag, std::size_t _ORD>
MATRICE_HOST_FINL auto _Invcomp_conv_impl<_Ty, _Tag, _ORD>::_Init(){
	const auto& _Pos = _Mybase::m_pos;
	const auto& _Ref = _Mybase::m_reference;
	auto[_L, _R, _U, _D] = conv_internal::_Square_range(_Pos, _Mybase::m_ksize>>1);

	_Mybase::_Myref = _Ref[0].block(_L, _R, _U, _D).eval();
	auto _Mean = _Mybase::_Myref.sum() / _Mybase::_Myref.size();
	_Mybase::_Myref = _Mybase::_Myref - _Mean;
	auto _SSD = sqrt((_Mybase::_Myref*_Mybase::_Myref).sum());
	if (_SSD < _Mybase::m_options)
		throw std::runtime_error("Bad value of variable _SSD in _Invcomp_conv_impl<...>::_Init().");
	_Mybase::_Myref = _Mybase::_Myref / _SSD;

	range<decltype(_L)> _Rx(_L, _R), _Ry(_U, _D);
	_Mybase::_Myjaco = conv_internal::_Op<1,6>::J(_Ref, _Pos, _Rx, _Ry);
	_Mybase::_Myhess = _Mybase::_Myjaco.t().mul(_Mybase::_Myjaco).reduce() * 2 / _SSD;

	_Mybase::_Myhess = _Mybase::_Mysolver.forward();
}
/**
 * \IC-GN solver kernel, includes computation of SDIs (Steepest Descent Images) and its dot-production with Error Image, backward substitution for solving dp and IC-update for warp parameters, which should be looped to convergence. 
 * \_Pars = {u, ux, uy, v, vx, vy} will be overwritten by the updated solution.
 */
template<typename _Ty, typename _Tag, std::size_t _ORD>
MATRICE_HOST_FINL auto _Invcomp_conv_impl<_Ty, _Tag, _ORD>::_Impl(param_type& _Pars) {
	// warp current image patch
	const auto _C = _Mybase::_Update_subset(_Pars);

	// error image expression
	auto _Diff_n = (_Mybase::_Myref - _Mybase::_Mycur)*_C;

	// computation of SDIs and its dot-production with error image
	stack_vector _Grad = (_Diff_n*_Mybase::_Myjaco).reduce().t();

	// solve $\Delta \mathbf{p}$
	_Grad = -1.*_Mybase::_Mysolver.backward(_Grad);

	// inverse-compositional update warp parameters $\mathbf{p}$
	_Pars = conv_internal::_Warp_update<_ORD>::inv(_Pars, _Grad);

	// update anchor position
	_Mybase::m_pos.x += _Pars[0], _Mybase::m_pos.y += _Pars[3];

	return std::make_tuple((_Diff_n*_Diff_n).sum(), (_Grad*_Grad).sum());
}

} MATRICE_ALGS_END