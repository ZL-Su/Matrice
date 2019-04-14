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
#include "../../core/matrix.h"
#include "../../core/vector.h"
#include "../../core/solver.h"
#include "../../private/math/_linear.h"
#include "../interpolation.h"
#include "_correlation_traits.h"

MATRICE_ALGS_BEGIN _DETAIL_BEGIN namespace corr {

struct _Correlation_options {
	size_t _Radius = 10;  //patch radius
	size_t _Maxits = 10;  //maximum iterations
	size_t _Stride =  7;  //node spacing
	float_t _Znssd= 0.4;  //correlation threshold

	template<typename _Ty, MATRICE_ENABLE_IF(is_floating_point_v<_Ty>)>
	static constexpr _Ty _Mytol = _Ty(1.0E-6); //iteration tolerance

	/**
	 * \retrieve range of the patch centered on point _Pos
	 */
	template<bool _Is_cutoff, typename _Ty, typename _Ret = conditional_t<_Is_cutoff, index_t, typename _Ty::value_t>>
	MATRICE_GLOBAL_INL auto range(const _Ty& _Pos) {
		using tuple_type = tuple<_Ret, _Ret, _Ret, _Ret>;
		if constexpr (_Is_cutoff) {
			const auto x = floor<_Ret>(_Pos.x), y = floor<_Ret>(_Pos.y);
			return tuple_type(x - _Radius, x + _Radius + 1, y - _Radius, y + _Radius + 1);
		}
		return tuple_type(_Pos.x-_Radius, _Pos.x+_Radius+1, _Pos.y - _Radius, _Pos.y + _Radius + 1);
	}

	template<typename _Pty, typename _Sty>
	MATRICE_HOST_INL auto range_check(const _Pty& _Pos, const _Sty& _Shape) {
		auto _TL = _Pos - _Radius; auto _RB = _Pos + _Radius;
		if (floor(_TL(0)) < 0 || floor(_TL(1)) < 0 ||
			_RB(0)>=get<0>(_Shape) || _RB(1)>=get<1>(_Shape))
			return false;
		return true;
	}
};

template<typename _Ty, typename _Itag, size_t _Order>
class _Corr_invcomp_optim {};
template<typename _Ty, typename _Itag, size_t _Order>
struct _Corr_solver_traits<_Corr_invcomp_optim<_Ty, _Itag, _Order>> {
	using value_type = _Ty;
	using itp_category = _Itag;
	static constexpr auto order = _Order;
};

template<typename _Derived> class _Corr_optim_base {
	using _Myt = _Corr_optim_base;
	using _Mydt = _Derived;
	using _Mytraits = _Corr_solver_traits<_Mydt>;
public:
	static constexpr auto DOF = _Mytraits::order*6;
	using value_type = typename _Mytraits::value_type;
	using matrix_type = Matrix<value_type>;
	using matrix_fixed = Matrix_<value_type, DOF, DOF>;
	using point_type = Vec2_<value_type>;
	using param_type = Vec_<value_type, DOF>;
	using option_type = _Correlation_options;
	using interp_type = typename interpolation<value_type, typename _Mytraits::itp_category>::type;
	using linear_solver = matrix_decomp<matrix_fixed, _TAG _Linear_spd_tag>;

	_Corr_optim_base(const interp_type& _Ref, const interp_type& _Cur, const point_type& _Pos, const option_type& _Opt) 
		:_Myref_itp(_Ref), _Mycur_itp(_Cur), _Mypos(_Pos), 
		 _Myopt(_Opt), _Mysolver(_Myhess) { 
		_Myref = this->_Init(); 
	}

	MATRICE_HOST_INL auto operator()(param_type& _p) {
		return static_cast<_Mydt*>(this)->_Update(_p);
	}

	// \for retrieve reference subset
	MATRICE_HOST_INL auto& f() const { return (_Myref); }

protected:
	MATRICE_HOST_INL auto& _Init();
	MATRICE_HOST_INL auto& _Warp(const param_type& _Pars);

	const interp_type& _Myref_itp;
	const interp_type& _Mycur_itp;

	option_type _Myopt;
	point_type _Mypos;
	value_type _Myissd;
	matrix_type  _Myref, _Mycur, _Myerr;
	matrix_type  _MyJaco;
	matrix_fixed _Myhess, _Hessian;
	linear_solver _Mysolver;
};

template<typename _Ty, typename _Itag>
class _Corr_invcomp_optim<_Ty, _Itag, 1> 
	: public _Corr_optim_base<_Corr_invcomp_optim<_Ty, _Itag, 1>>
{
	using _Myt = _Corr_invcomp_optim;
	using _Mybase = _Corr_optim_base<_Myt>;
public:
	static constexpr auto order = compile_time_size<>::val_1;
	using typename _Mybase::option_type;
	using typename _Mybase::param_type;
	using typename _Mybase::point_type;
	using typename _Mybase::value_type;

	_Corr_invcomp_optim(const _Myt& _Other) = delete;
	_Corr_invcomp_optim(_Myt&& _Other) = delete;
	template<typename... _Args>
	_Corr_invcomp_optim(const _Args&..._args)
		: _Mybase(_args...) {}

	/**
	 *\brief Solve normal equations to update parameters
	 *\param [_P] warp parameter vector
	 */
	MATRICE_HOST_INL auto _Update(param_type& _P);

	/**
	 *\brief Compute steepest descent image terms
	 */
	MATRICE_HOST_INL auto& _Diff();
};

_DETAIL_END } MATRICE_ALGS_END

DGE_MATRICE_BEGIN
struct xcorr_optim {
	using option = algs::detail::corr::_Correlation_options;
	template<typename _Ty, typename _Itag, size_t _Order=1>
	using ic = algs::detail::corr::_Corr_invcomp_optim<_Ty, _Itag, _Order>;
};
DGE_MATRICE_END

#include "inline\_optim.inl"