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
#include "../_icgn_impl.h"
#include "../../../private/_range.h"

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
 * \estimate Jacobian: df/d\mathbf{p} = $\partial f / \partial p_i$
 */
template<> struct _Op<1, 6> {
	template<typename _Mulmty, typename _Poty, typename _Rgty>
	static auto J(const _Mulmty& _Ref, const _Poty& _Pos, const _Rgty& _Rx, const _Rgty& _Ry) {
		using value_type = typename _Mulmty::value_type;

		auto [_Dfdx, _Dfdy] = _Ref.view_n<2,1>(_Rx, _Ry);
		const auto _Size = max(_Rx.size(), _Ry.size());

		tensor<value_type, 1, 6> _Ret(_Size, _Size);
		for (auto y : _Ry) { auto _Dy = y - _Pos.y;
			const auto _y = y - _Ry.begin();
			for (auto x : _Rx) { auto _Dx = x - _Pos.x;
				const auto _x = x - _Rx.begin();
				const auto &_Fx = _Dfdx[y][x], &_Fy = _Dfdy[y][x];
				_Ret(_y, _x) = { _Fx, _Fy, _Fx*_Dx, _Fx*_Dy, _Fy*_Dx, _Fy*_Dy };
			}
		}
		return std::forward<decltype(_Ret)>(_Ret);
	};
};
}

/**
 * \update deformed subset for each interation
 */
template<typename _Derived> MATRICE_HOST_FINL 
auto _Iterative_conv_base<_Derived>::_Update_subset(const param_type& _P) {
#define _WITHIN_RANGE_OF_REFIMG \
	if (_Bu>=0&&_Bl>=0 && _Bd<m_reference[0].rows()&&_Br<m_reference[0].cols())

	using _Interp_t = conv_internal::_Op<myt_traits_t::interp, 1, decltype(m_coeff), stack_vector>;

	const auto _Radius = static_cast<int>(m_options());
	const auto& _Center = m_pos;
	const auto[_L, _R, _U, _D] = conv_internal::_Square_range(_Center, _Radius);

	stack_vector _Buf_x{ 1. }, _Buf_y{ 1. };
	_Interp_t _Op(m_coeff, _Buf_x, _Buf_y);

	auto _Mean = zero<value_type>::value;
	for (auto y = _U; y < _D; ++y) {
		auto _Dy = static_cast<value_type>(y - _Center.y);
		auto _Y = y + _P[3] + _P[5] * _Dy;
		auto _X = _P[0] + _P[2] * _Dy;
		for (auto x = _L; x < _R; ++x) {
			auto _Dx = static_cast<value_type>(x - _Center.x);
			_X += x + _P[1] * _Dx; _Y += _P[4] * _Dx;

			auto _Ix = floor<int>(_X), _Iy = floor<int>(_Y);
			auto _Bu = _Iy - 2, _Bd = _Iy + 3;
			auto _Bl = _Ix - 2, _Br = _Ix + 3;
			_WITHIN_RANGE_OF_REFIMG{
				auto _Delta_x = _X - static_cast<value_type>(_Ix);
				auto _Delta_y = _Y - static_cast<value_type>(_Iy);

				conv_internal::_Delta_pow_n(_Buf_x, _Delta_x);
				conv_internal::_Delta_pow_n(_Buf_y, _Delta_y);

				_Ix = _Bl * param_type::Size;
				_Iy = _Bu * param_type::Size;
				_Mean += m_current[y - _U][x - _L] = _Op(_Ix, _Iy);
			}
			else m_current[y - _U][x - _L] = zero<value_type>::value;
		}
	}
	_Mean /= static_cast<value_type>(multiply(_R - _L, _D - _U));

	auto _Diff_exp = m_current - _Mean;
	auto _SSD = sqrt((_Diff_exp*_Diff_exp).sum());

	return std::make_tuple(range(_L, _R), range(_U, _D), _Mean, _SSD);
}

// \IC-GN solver: _Pars={u, ux, vx, v, uy, vy} will be overwritten by the updated solution
template<typename _Ty, std::size_t _Intp, std::size_t _Order> MATRICE_HOST_FINL 
auto _Invcomp_conv_impl<_Ty, _Intp, _Order>::_Solver_impl(param_type& _Pars) {

	// \step 1: calc. Jacobian and Hessian;
	// \step 2: iteratively update warp pars.

	const auto& _Ref  = _Mybase::m_reference;

	const auto[_Rx, _Ry, _G_mean, _G_ssd] = _Mybase::_Update_subset(_Pars);
	if (_G_ssd < _Mybase::m_options) throw std::runtime_error("Bad value of variable _G_ssd in _Invcomp_conv_impl<_Ty, _Intp, _Order>::_Solver_impl(...).");
	
	// \estimate df/dp
	const auto _Jaco = conv_internal::_Op<_Order, _Mybase::DOF>::J(_Ref, _Mybase::m_pos, _Rx, _Ry);

	auto _Hess = _Jaco.t().mul(_Jaco).eval();

	auto _Diff_f_exp = _Ref[0].block(_Rx.begin(), _Rx.end(), _Ry.begin(), _Ry.end()).eval() - m_favg;
	auto _Diff_g_exp = m_current - _G_mean;
	auto _Ndiff_exp = _Diff_f_exp * (1 / m_fssd) - _Diff_g_exp * (1 / _G_ssd);

	auto _Coef = zero<value_type>::value;
	/*for (auto y : _Ry) {
		for (auto x : _Rx) {
			auto _Diff_f = _Diff_f_exp(x, y);
			_Diff_f /= m_fssd;
			auto _Diff_g = _Diff_g_exp(x - _Rx.front(), y - _Ry.front());
			_Diff_g /= _G_ssd;

			auto _Diff = _Diff_f - _Diff_g;

		}
	}*/

	stack_vector _Grad = (_Ndiff_exp*_Jaco).sum()*(2/_G_ssd);

#undef _WITHIN_RANGE_OF_REFIMG
}


} MATRICE_ALGS_END