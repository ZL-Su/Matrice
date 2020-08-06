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
#include <variant>
#include "core/solver.h"
#include "private/math/_linear.h"
#include "_correlation_traits.h"
#include "../interpolation.h"

MATRICE_ALG_BEGIN(corr)
_DETAIL_BEGIN

struct _Correlation_options {
	size_t _Radius = 15;  //patch radius
	size_t _Maxits = 50;  //maximum iterations
	size_t _Stride =  7;  //node spacing
	float_t _Znssd = 0.6; //correlation threshold
	float_t _Coeff = 0.7; //damping coefficient

	template<typename _Ty = default_type, MATRICE_ENABLE_IF(is_floating_point_v<_Ty>)>
	static constexpr _Ty _Mytol = _Ty(1.0E-6); //iteration tolerance

	/**
	 * \retrieve range of the patch centered on point _Pos
	 */
	template<bool _Is_cutoff, typename _Ty, 
		typename _Ret = conditional_t<_Is_cutoff, index_t, typename _Ty::value_t>>
	MATRICE_GLOBAL_INL auto range(const _Ty& _Pos)->tuple<_Ret, _Ret, _Ret, _Ret> {
		using tuple_type = tuple<_Ret, _Ret, _Ret, _Ret>;
		if constexpr (_Is_cutoff) {
			const auto x = floor<_Ret>(_Pos.x), y = floor<_Ret>(_Pos.y);
			return tuple_type(x - _Radius, x + _Radius + 1, y - _Radius, y + _Radius + 1);
		}
		return tuple_type(_Pos.x-_Radius, _Pos.x+_Radius+1, _Pos.y - _Radius, _Pos.y + _Radius + 1);
	}

	template<typename _Pty, typename _Sty>
	MATRICE_HOST_INL bool range_check(const _Pty& _Pos, const _Sty& _Shape) {
		auto _TL = _Pos - _Radius; auto _RB = _Pos + _Radius;
		if (floor(_TL(0)) < 0 || floor(_TL(1)) < 0 ||
			_RB(0)>=get<0>(_Shape) || _RB(1)>=get<1>(_Shape))
			return false;
		return true;
	}
};

// \solver tag declarations
template<size_t _Order> struct _Alg_icgn : _TAG _Solver_tag {
	static constexpr auto order = _Order;
	static_assert(order < 3, "_Order must be 0, 1 or 2");
}; // inverse-compositional GN
template<size_t _Order> struct _Alg_fagn : _TAG _Solver_tag {
	static constexpr auto order = _Order;
	static_assert(order < 3, "_Order must be 0, 1 or 2");
}; // forward-additional GN

// \forward solver declaration 
template<typename _Ty, typename _Itag, typename _Atag>
class _Corr_solver_impl {};

// \update on warp parameters
template<typename _Alg_tag> struct _Param_update_strategy {};
template<typename _Itp_tag> struct _Corr_border_size {};

template<typename _Ty, typename _Itag, typename _Atag>
struct _Corr_solver_traits<_Corr_solver_impl<_Ty, _Itag, _Atag>> {
	using value_type = _Ty;
	using itp_category = _Itag;
	using alg_category = _Atag;
	using interpolator = interpolation<value_type, itp_category>;
	using update_strategy = _Param_update_strategy<alg_category>;
	using border_size = _Corr_border_size<itp_category>;
	static constexpr auto order = alg_category::order;
};

// \correlation optimizer base class
template<typename _Derived> class _Corr_optim_base {
	using _Myt = _Corr_optim_base;
	using _Mydt = _Derived;
protected:
	using _Mytraits = _Corr_solver_traits<_Mydt>;
public:
	static constexpr auto npar = conditional_size_v<
		_Mytraits::order == 0, 2, _Mytraits::order * 6>;
	using value_type = typename _Mytraits::value_type;
	using matrix_type = Matrix<value_type>;
	using matrix_fixed = Matrix_<value_type, npar, npar>;
	using point_type = Vec2_<value_type>;
	using param_type = Vec_<value_type, npar>;
	using jacob_type = Matrix_<value_type, ::dynamic, npar>;
	using vector_type = Matrix_<value_type, ::dynamic, 1>;
	using options_type = _Correlation_options;
	using interp_type = typename _Mytraits::interpolator;
	using update_strategy = typename _Mytraits::update_strategy;
	using linear_solver = matrix_decomp<matrix_fixed, _TAG _Linear_spd_tag>;
	using rect_type = rect<size_t>;

	/**
	 *\brief This module is image-wise thread-safe. It allows us to eval. relative deformation of each point between current and reference state.
	 *\param _Ref reference interpolated image;
	 *\param _Cur current interpolated image;
	 *\param _Opt options for the optimizer.
	 */
	_Corr_optim_base(
		const interp_type& _Ref,
		const interp_type& _Cur,
		const options_type& _Opt
	) : _Myopt(_Opt),
		_Myref_ptr(std::make_shared<interp_type>(_Ref)),
		_Mycur_ptr(std::make_shared<interp_type>(_Cur)),
		_Mysolver(_Myhess), _Mysize(_Opt._Radius*2+1),
		_Myjaco(sq(_Mysize)), _Mydiff(sq(_Mysize)){
	}

	/**
	 *\brief set a reference point being estimated. 
	 *\param [_Pos] a reference position for parameter estimation;
	 */
	template<typename _Ty>
	MATRICE_HOST_INL _Myt& init(const _Ty& _x, const _Ty& _y) {
		_Mypos.x = _x, _Mypos.y = _y;
		/*_Myref = */this->_Cond();
		return (*this);
	}
	MATRICE_HOST_INL _Myt& init(const point_type& _ref_pos) {
		_Mypos.x = _ref_pos.x, _Mypos.y = _ref_pos.y;
		/*_Myref = */this->_Cond();
		return (*this);
	}

	/**
	 *\brief Integer pixel level search.
	 *\param roi region of interest
	 */
	MATRICE_HOST_INL point_type guess(rect_type&& roi) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(!_Myref.empty, "Call init(...) before calling method guess(...).");
#endif
		return this->_Guess(roi);
	}

	/**
	 *\brief solve Gauss-Newton normal equations
	 *\param (_pars) warp parameters.
	 *\return [znssd, error := dot(dp, dp)]
	 */
	MATRICE_HOST_INL auto operator()(param_type& _pars) {
		return (this->_Solve(_pars));
	}

	// \for retrieving reference subset
	MATRICE_HOST_INL decltype(auto)get_refpatch() const noexcept { 
		return (_Myref); 
	}
	MATRICE_HOST_INL constexpr decltype(auto)options()const noexcept{
		return (_Myopt);
	}
	MATRICE_HOST_INL constexpr decltype(auto)options() noexcept {
		return (_Myopt);
	}
	MATRICE_HOST_INL decltype(auto)refpos() noexcept { 
		return _Mypos; 
	}
	MATRICE_HOST_INL decltype(auto)refpos() const noexcept { 
		return _Mypos; 
	}

protected:
	///<methods>

	/**
	 *\brief Build refpatch, and buffs to hold curpatch and diffs.
	 *\return the refpatch: this->_Myref.
	 */
	MATRICE_HOST_INL auto _Cond()->matrix_type&;

	/**
	 *\brief Build refpatch, and buffs to hold curpatch and diffs.
	 *\return the refpatch: this->_Myref.
	 */
	MATRICE_HOST_INL auto _Guess(rect_type roi)->point_type;

	/**
	 *\brief Solve new parameters
	 *\param [_Par] in: old parameters, output: new parameters
	 */
	MATRICE_HOST_INL auto _Solve(param_type& _Pars);

	///</methods>

	///<fields>
	diff_t       _Mysize;
	options_type  _Myopt;
	point_type   _Mypos;
	matrix_type  _Myref, _Mycur;
	vector_type  _Mydiff;
	jacob_type   _Myjaco;
	matrix_fixed _Myhess;

	linear_solver _Mysolver;

	// \hold interpolated reference image
	std::shared_ptr<interp_type> _Myref_ptr;
	// \hold interpolated current image 
	std::shared_ptr<interp_type> _Mycur_ptr;
	///</fields>
};

// \zero-order inverse compositional gauss-newton algrithom
template<typename _Ty, typename _Itag>
class _Corr_solver_impl<_Ty, _Itag, _Alg_icgn<0>>
	: public _Corr_optim_base<_Corr_solver_impl<_Ty, _Itag, _Alg_icgn<0>>>
{
	using _Myt = _Corr_solver_impl;
	using _Mybase = _Corr_optim_base<_Myt>;
	using typename _Mybase::_Mytraits;
public:
	static constexpr auto order = _Mytraits::order;
	using typename _Mybase::options_type;
	using typename _Mybase::param_type;
	using typename _Mybase::point_type;
	using typename _Mybase::value_type;

	_Corr_solver_impl(const _Myt& _Other) = delete;
	_Corr_solver_impl(_Myt&& _Other) = delete;
	template<typename... _Args>
	_Corr_solver_impl(const _Args&..._args) : _Mybase(_args...) {}

	/**
	 *\brief Evaluate Jacobian contribution at ref. patch.
	 */
	MATRICE_HOST_INL auto& _Diff();

	/**
	 *\brief Construct current image patch.
	 *\param [_Par] displacement {u, v}.
	 */
	MATRICE_HOST_INL auto& _Warp(const param_type& _Par);

private:
	using _Mybase::_Myref_ptr;
	using _Mybase::_Mycur_ptr;
	using _Mybase::_Mycur;
	using _Mybase::_Myopt;
	using _Mybase::_Mypos;
};

// \1-st order inverse compositional gauss-newton algrithom
template<typename _Ty, typename _Itag>
class _Corr_solver_impl<_Ty, _Itag, _Alg_icgn<1>>
	: public _Corr_optim_base<_Corr_solver_impl<_Ty, _Itag, _Alg_icgn<1>>>
{
	using _Myt = _Corr_solver_impl;
	using _Mybase = _Corr_optim_base<_Myt>;
	using typename _Mybase::_Mytraits;
public:
	static constexpr auto order = _Alg_icgn<1>::order;
	using typename _Mybase::options_type;
	using typename _Mybase::param_type;
	using typename _Mybase::point_type;
	using typename _Mybase::value_type;

	_Corr_solver_impl(const _Myt& _Other) = delete;
	_Corr_solver_impl(_Myt&& _Other) = delete;
	template<typename... _Args>
	_Corr_solver_impl(const _Args&..._args)
		: _Mybase(_args...) {}

	/**
	 *\brief evaluate Jacobian contribution at refpatch.
	 */
	MATRICE_HOST_INL auto& _Diff();
	/**
	 *\brief Construct current image patch.
	 *\param [_Par] {u, dudx, dudy, v, dvdx, dvdy}
	 */
	MATRICE_HOST_INL auto& _Warp(const param_type& _Par);

private:
	using _Mybase::_Myref_ptr;
	using _Mybase::_Mycur_ptr;
	using _Mybase::_Mycur;
	using _Mybase::_Myopt;
	using _Mybase::_Mypos;
};

_DETAIL_END
MATRICE_ALG_END(corr)

#include "inline/_optim.inl"