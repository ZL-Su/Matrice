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
***********************************************************************/
#pragma once

#include "../_optim.h"
#ifdef MATRICE_SIMD_ARCH
#include "arch/ixpacket.h"
#endif
#include "../../similarity.h"
#include "private/math/_linear_solver.hpp"

MATRICE_ALG_BEGIN(corr)
_DETAIL_BEGIN

// \retrieve border size for each interpolation alg.
template<> struct _Corr_border_size<bilerp_tag> {
	static constexpr auto lower = 1, upper = 1;
};
template<> struct _Corr_border_size<bicerp_tag> {
	static constexpr auto lower = 1, upper = 2;
};
template<> struct _Corr_border_size<biqerp_tag> {
	static constexpr auto lower = 2, upper = 3;
};
template<> struct _Corr_border_size<biserp_tag> {
	static constexpr auto lower = 3, upper = 4;
};

// \specializations of parameter update strategy.
template<> struct _Param_update_strategy<_Alg_icgn<0>> {
	template<typename _Ty>
	static MATRICE_GLOBAL_FINL _Ty& eval(_Ty& x, const _Ty& y) {
		x[0] -= y[0], x[1] -= y[1];
		return (x);
	}
};
template<> struct _Param_update_strategy<_Alg_icgn<1>> {
	template<typename _Ty>
	static MATRICE_GLOBAL_FINL _Ty& eval(_Ty& x, const _Ty& y) {
		constexpr auto _One = one<typename _Ty::value_t>;

		const auto y1p1 = _One + y[1], y5p1 = _One + y[5];
		const auto _Inv_of_det = _One / (y1p1*y5p1 - y[2]*y[4]);
		const auto a = y[2] * y[3] - y[0] * y5p1;
		const auto b = y[0] * y[4] - y[3] * y1p1;

		const auto x1p1 = _One + x[1], x5p1 = _One + x[5];
		const auto u  = _Inv_of_det*(x1p1*a    + x[2]*b) + x[0];
		const auto ux = _Inv_of_det*(x1p1*y5p1 - x[2]*y[4]) - _One;
		const auto uy = _Inv_of_det*(x[2]*y1p1 - y[2]*x1p1);
		const auto v  = _Inv_of_det*(x[4]*a    + x5p1*b) + x[3];
		const auto vx = _Inv_of_det*(x[4]*y5p1 - y[4]*x5p1);
		const auto vy = _Inv_of_det*(x5p1*y1p1 + x[4]*y[2]) - _One;

		x[0] = u, x[1] = ux, x[2] = uy, x[3] = v, x[4] = vx, x[5] = vy;

		return (x);
	}
};

#pragma region <-- base class implementation -->
template<typename _Derived> MATRICE_HOST_INL 
auto _Corr_optim_base<_Derived>::_Cond()->matrix_type& {
	// \create buf.s to hold reference and current patchs
	_Myref.create(_Mysize, _Mysize);
	_Mycur.create(_Mysize, _Mysize);

	// \fill reference image patch.
	const auto& _Data = _Myimref->data();
	const auto[_Rows, _Cols] = _Data.shape().tile();
	const auto[_L, _R, _U, _D] = _Myopt.range<true>(_Mypos);
	auto _Mean = zero<value_type>;
	if (_Mypos.x == floor(_Mypos.x) && _Mypos.y == floor(_Mypos.y)) {
		for (auto y = _U; y < _D; ++y) {
			if (y >= 0 && y < _Rows) {
				const auto _Dp = _Data[y]; 
				auto _Rp = _Myref[y - _U];
				for (auto x = _L; x < _R; ++x) {
					if (x >= 0 && x < _Cols)
						_Mean += _Rp[x - _L] = _Dp[x];
				}
			}
		}
	} else { //fill _Myref with interpolation for non-integer targets.
		for (auto y = _U; y < _D; ++y) {
			if (y >= 0 && y < _Rows) {
				auto _Rp = _Myref[y - _U];
				for (auto x = _L; x < _R; ++x) {
					if (x >= 0 && x < _Cols)
						_Mean += _Rp[x - _L] = (*_Myimref)(x, y);
				}
			}
		}
	}

	// \comp. Jacobian and Hessian
	static_cast<_Derived*>(this)->_Diff();
#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
	_Myhess = _Myjaco.t().mul_inplace(_Myjaco);
#else
	_Myhess = _Myjaco.t().mul(_Myjaco);
#endif

	// \zero mean normalization.
	_Myref = _Myref - (_Mean /= _Mysize*_Mysize);
	const auto _Issd = one<value_type> / sqrt(sq(_Myref).sum());
	_Myref = _Myref * _Issd;
	
	_Myhess = _Myhess * _Issd;
	_Myhess = _Myhess + matrix_fixed::diag(_Myopt._Coeff);
	internal::_Lak_adapter<spt>(_Myhess.view());
	
	// ad hoc for robust eval.
	_Myissd = _Issd;

	return (_Myref);
}

template<typename _Derived> MATRICE_HOST_INL 
auto _Corr_optim_base<_Derived>::_Guess(rect_type roi)->point_type {
	const auto _Start = roi.begin(), _End = roi.end();
	const auto _Data = _Myimcur->data();
	const auto _Stride = rect_type::value_type(_Myopt._Radius>>1);

	//narrow the ROI if it hits the boundaries of the image
	const auto _Off{ _Myopt._Radius + 1 };
	transforms::clamp<rect_type::value_type> _Clamp(_Off, _Data.cols()-_Off);
	_Start[0] = _Clamp(_Start.x), _End[0] = _Clamp(_End.x);
	_Clamp._Myupper = _Data.rows() - _Off;
	_Start[1] = _Clamp(_Start.y), _End[1] = _Clamp(_End.y);

	zncc_metric_t<value_type> zncc(_Myref.data(), _Myopt._Radius);

	point_type _Pos;
	value_type _Max{ 0 };
	for (auto y = _Start.y; y < _End.y; y += _Stride) {
		for (auto x = _Start.x; x < _End.x; x += _Stride) {
			_Mycur = _Data.block(x, y, _Myopt._Radius);
			const auto _Coeff  = zncc.eval(_Mycur.data());
			if (_Coeff > _Max) {
				_Max = _Coeff;
				_Pos.x = x, _Pos.y = y;
			}
		}
	}
	const auto _Hs = _Stride >> 1;
	for (auto y = size_t(_Pos.y) - _Hs; y < size_t(_Pos.y) + _Hs; ++y){
		for (auto x = size_t(_Pos.x) - _Hs; x < size_t(_Pos.x) + _Hs; ++x){
			_Mycur = _Data.block(x, y, _Myopt._Radius);
			const auto _Coeff = zncc.eval(_Mycur.data());
			if (_Coeff > _Max) {
				_Max = _Coeff;
				_Pos.x = x, _Pos.y = y;
			}
		}
	}
	return _Pos;
}

template<typename _Derived> MATRICE_HOST_INL
auto _Corr_optim_base<_Derived>::_Solve(param_type& Par) {
	// \warp current image patch.
	_Mycur = static_cast<_Derived*>(this)->_Warp(Par);

	// \error map
	_Mydiff = _Mycur - _Myref;
#ifdef MATRICE_DEBUG
	matrix_type _Error_map(_Mycur.shape(), _Mydiff.data());
#endif

	// \steepest descent parameter update
#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
	param_type _Sdp = _Myjaco.t().mul_inplace(_Mydiff);
#else
	param_type _Sdp = _Myjaco.t().mul(_Mydiff);
#endif

	// \solve update to the warp parameter vector.
	internal::_Bwd_adapter<spt>(_Myhess.view(), _Sdp.view());

	// \inverse composition to update param.
	Par = update_strategy::eval(Par, _Sdp);

	// \report least square correlation coeff. and param. error.
	return std::make_tuple(sq(_Mydiff).sum(), _Sdp.dot(_Sdp));
}
template<typename _Derived> MATRICE_HOST_INL 
auto _Corr_optim_base<_Derived>::robust_sol(param_type& Par)
{
	// \warp current image patch.
	_Mycur = static_cast<_Derived*>(this)->_Warp(Par);

	// \error map
	_Mydiff = _Mycur - _Myref;
#ifdef MATRICE_DEBUG
	const matrix_type _Error_map(_Mycur.shape(), _Mydiff.data());
#endif

	// \compute robustness and weight Jacobian
	for (auto i = 0; i < _Myweight.size(); ++i) {
		const auto wi = _Myloss.phi(sq(_Mydiff(i)));
		_Myweight(i) = wi;
		for (auto j = 0; j < _Myjaco.cols(); ++j) {
			_Myjaco_tw.cview(i)(j) = _Myjaco[i][j] * wi;
		}
	}

	// \compute weighted Gauss-Newton Hessian matrix
	_Myhess = _Myjaco_tw.mul(_Myjaco);
	_Myhess = _Myhess * _Myissd;

	// \steepest descent parameter update
	param_type _Sdp = _Myjaco_tw.mul(_Mydiff);

	// \solve update to the warp parameter vector.
	const auto solver = make_linear_solver<lud>(_Myhess);
	const auto _Dp = solver.solve(_Sdp);

	// \inverse composition to update param.
	Par = update_strategy::eval(Par, _Dp);

	auto _Loss = value_type(0);
	for (auto i = 0; i < _Mydiff.size(); ++i) {
		_Loss += _Myloss.rho(sq(_Mydiff(i)));
	}

	// \report least square correlation coeff., param. error, and loss.
	return std::make_tuple(sq(_Mydiff).sum(), _Loss, _Dp.dot(_Dp));
}
#pragma endregion

#pragma region <-- derived class implementation -->
#define _IF_WITHIN_RANGE(_OP) \
	if (y-_Mytraits::border_size::lower>=0 && \
		 x-_Mytraits::border_size::lower>=0 && \
		 y+_Mytraits::border_size::upper<_Rows && \
		 x+_Mytraits::border_size::upper<_Cols) \
		_OP; \
   else _Mycur(r,c) = zero<value_type>;

// \specialization of warp. to eval. current patch.
template<typename _Ty, typename _Itag> MATRICE_HOST_INL
auto& _Corr_solver_impl<_Ty, _Itag, _Alg_icgn<0>>::_Warp(const param_type& _Par) {

	const auto[_Rows, _Cols] = (*_Myimcur)().shape().tile();
	const auto _Radius = static_cast<index_t>(_Myopt._Radius);

	const auto _Cur_x = _Mypos[0] + _Par[0];
	const auto _Cur_y = _Mypos[1] + _Par[1];

	auto _Mean = zero<value_type>;
	for (diff_t j = -_Radius, r = 0; j <= _Radius; ++j, ++r) {
		const auto y = _Cur_y + static_cast<value_type>(j);
		auto _Rc = _Mycur[r];
		for (diff_t i = -_Radius, c = 0; i <= _Radius; ++i, ++c) {
			const auto x = static_cast<value_type>(i) + _Cur_x;
			_IF_WITHIN_RANGE(_Mean += _Rc[c] = (*_Myimcur)(x, y));
		}
	}

	_Mycur = _Mycur - (_Mean /= sq(_Mybase::_Mysize));
	const auto _Issd = one<value_type>/sqrt(sq(_Mycur).sum());
	_Mycur = _Mycur * _Issd;
	return (_Mycur);
}
template<typename _Ty, typename _Itag> MATRICE_HOST_INL
auto& _Corr_solver_impl<_Ty, _Itag, _Alg_icgn<1>>::_Warp(const param_type& _Par) {
	const auto[_Rows, _Cols] = _Myimcur->data().shape().tile();
	const auto _Radius = static_cast<index_t>(_Myopt._Radius);

	const auto u = _Par[0], dudx = _Par[1], dudy = _Par[2];
	const auto v = _Par[3], dvdx = _Par[4], dvdy = _Par[5];
	const auto _Cur_x = _Mypos[0] + u, _Cur_y = _Mypos[1] + v;

	auto _Mean = zero<value_type>;
	for (diff_t j = -_Radius, r = 0; j <= _Radius; ++j, ++r) {
		const auto dy = static_cast<value_type>(j);
		const auto tx = _Cur_x + dudy * dy;
		const auto ty = _Cur_y + (1 + dvdy) * dy;
		auto _Rc = _Mycur[r];
		for (diff_t i = -_Radius, c = 0; i <= _Radius; ++i, ++c) {
			const auto dx = static_cast<value_type>(i);
			const auto x = (1 + dudx)* dx + tx;
			const auto y = dvdx * dx + ty;
			_IF_WITHIN_RANGE(_Mean += _Rc[c] = (*_Myimcur)(x, y));
		}
	}

	_Mycur = _Mycur - (_Mean /= sq(_Mybase::_Mysize));
	const auto _Issd = one<value_type> / sqrt(sq(_Mycur).sum());
	_Mycur = _Mycur * _Issd;

	return (_Mycur);
}

// \specialization of diff. to eval. Jacobian contributions.
template<typename _Ty, typename _Itag> MATRICE_HOST_INL
auto& _Corr_solver_impl<_Ty, _Itag, _Alg_icgn<0>>::_Diff() {
	const auto _Size = _Mybase::_Mysize;
	const auto[_L, _R, _U, _D] = _Myopt.range<true>(_Mypos);
	const auto _Off = -static_cast<value_type>(_Myopt._Radius);

	for (index_t iy = _U, j = 0; iy < _D; ++iy, ++j) {
		const auto y = _Mypos.y + _Off + j;
		for (index_t ix = _L, i = 0; ix < _R; ++ix, ++i) {
			const auto x = _Mypos.x + _Off + i;
			const auto[dfdx, dfdy] = _Myimref->grad({ x, y });

			auto q = _Mybase::_Myjaco[j * _Size + i];
			q[0] = dfdx, q[1] = dfdy;
		}
	}

	return (_Mybase::_Myjaco);
}
template<typename _Ty, typename _Itag> MATRICE_HOST_INL
auto& _Corr_solver_impl<_Ty, _Itag, _Alg_icgn<1>>::_Diff() {
	const auto _Size = _Mybase::_Mysize;
	const auto[_L, _R, _U, _D] = _Myopt.range<true>(_Mypos);
	const auto _Off = -static_cast<value_type>(_Myopt._Radius);

	for (index_t iy = _U, j = 0; iy < _D; ++iy, ++j) {
		const auto dy = _Off + j;
		const auto y = _Mypos.y + dy;
		for (index_t ix = _L, i = 0; ix < _R; ++ix, ++i) {
			const auto dx = _Off + i;
			const auto x = _Mypos.x + dx;
			const auto[dfdx, dfdy] = _Myimref->grad({ x, y });

			auto q = _Mybase::_Myjaco[j * _Size + i];
			q[0] = dfdx, q[1] = dfdx * dx, q[2] = dfdx * dy;
			q[3] = dfdy, q[4] = dfdy * dx, q[5] = dfdy * dy;
		}
	}

	return (_Mybase::_Myjaco);
}

template<typename _Ty, typename _Itag> MATRICE_HOST_INL 
auto _Corr_solver_impl<_Ty, _Itag, _Alg_icgn<1>>::_Diff(value_type cx, value_type cy)
{
	const size_t _Size = _Mybase::_Mysize;
	const auto [_L, _R, _U, _D] = _Myopt.range<true>(point_type{ cx,cy });
	const auto _Off = -static_cast<value_type>(_Myopt._Radius);

	typename _Mybase::jacob_type _Jacobian(sq(_Size));
	for (index_t iy = _U, j = 0; iy < _D; ++iy, ++j) {
		const auto dy = _Off + j;
		const auto y = cy + dy;
		for (index_t ix = _L, i = 0; ix < _R; ++ix, ++i) {
			const auto dx = _Off + i;
			const auto x = cx + dx;
			const auto [dfdx, dfdy] = _Myimref->grad({ x, y });

			auto q = _Jacobian[j * _Size + i];
			q[0] = dfdx, q[1] = dfdx * dx, q[2] = dfdx * dy;
			q[3] = dfdy, q[4] = dfdy * dx, q[5] = dfdy * dy;
		}
	}

	return move(_Jacobian);
}

#undef _IF_WITHIN_RANGE
#pragma endregion

_DETAIL_END

/// <summary>
/// \brief Evaluate performance in terms of mean and standard deviation errors.
/// </summary>
/// <typeparam name="_FwdIt">Require forward iterator type</typeparam>
/// <param name="_Begin">: Start iterator of measurements.</param>
/// <param name="_End'">: End iterator of measurements.</param>
/// <param name="_Gt">: Ground-truth.</param>
/// <returns>Tuple with {mean, SD}</returns>
template<typename _FwdIt>
MATRICE_GLOBAL_INL auto eval_perf(_FwdIt _Begin, _FwdIt _End, default_type _Gt) noexcept {
	decltype(_Gt) _Mean = 0;
	for (auto _It = _Begin; _It != _End; ++_It) {
		_Mean += *_It;
	}
	const auto _Count = _End - _Begin;
	_Mean = _Mean / _Count - _Gt;

	decltype(_Mean) _Stdv = 0;
	for (auto _It = _Begin; _It != _End; ++_It) {
		_Stdv += sq(abs(*_It - _Gt) - _Mean);
	}
	_Stdv = sqrt(_Stdv / (_Count - 1));

	return tuple{ _Mean, _Stdv };
}

MATRICE_ALG_END(corr)