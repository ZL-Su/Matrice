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
#include <variant>
#include "_correlation_traits.h"
#include "../interpolation.h"

MATRICE_ALG_BEGIN(corr)
_DETAIL_BEGIN

// \solver tag declarations
/// <summary>
/// \brief IC-GN algorithm tag for internal solver
/// </summary>
template<size_t _Order> 
struct _Alg_icgn : _TAG _Solver_tag {
	static constexpr auto order = _Order;
	static_assert(order < 3, "_Order must be 0, 1 or 2");
};

/// <summary>
/// \brief FA-NR algorithm tag for internal solver
/// </summary>
template<size_t _Order> 
struct _Alg_fanr : _TAG _Solver_tag {
	static constexpr auto order = _Order;
	static_assert(order < 3, "_Order must be 0, 1 or 2");
};

// \forward solver declaration 
template<typename _Ty, typename _Itag, typename _Atag>
class _Corr_solver_impl {};

// \update on warp parameters
template<typename _Alg_tag> struct _Param_update_strategy {};
template<typename _Itp_tag> struct _Corr_border_size {};

/// <summary>
/// \brief TRAITS CLASS, correlation optimizer traits
/// </summary>
/// <typeparam name="_Ty">Scalar, primitive value type</typeparam>
/// <typeparam name="_Itag">Interpolation function tag</typeparam>
/// <typeparam name="_Atag">Internal solver tag</typeparam>
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

/// <summary>
/// \brief CLASS, correlation solver options
/// </summary>
struct _Correlation_options {
	size_t _Radius = 15;  //patch radius
	size_t _Maxits = 50;  //maximum iterations
	size_t _Stride = 7;  //node spacing
	float_t _Znssd = 0.6; //correlation threshold
	float_t _Coeff = 0.0; //damping coefficient
	mutable float_t _Mytol = 1.0e-6;

	/**
	 * \brief Setter and geter of subset radius.
	 */
	decltype(auto) radius() const noexcept {
		return (_Radius);
	}
	decltype(auto) radius() noexcept {
		return (_Radius);
	}

	/**
	 * \brief Setter and geter of max iterations.
	 */
	decltype(auto) maxiters() const noexcept {
		return (_Maxits);
	}
	decltype(auto) maxiters() noexcept {
		return (_Maxits);
	}

	/**
	 * \brief Setter and geter of correlation threshold.
	 */
	decltype(auto) threshold() const noexcept {
		return (_Znssd);
	}
	decltype(auto) threshold() noexcept {
		return (_Znssd);
	}

	/**
	 * \brief Geter of correlation threshold with base 1E-6.
	 */
	template<typename _Ty = default_type>
	decltype(auto) tol(_Ty const scale) const noexcept {
		return (scale * _Mytol);
	}

	/**
	 * \brief Get range of the patch centered at point '_Pos'.
	 */
	template<bool _Is_cutoff, typename _Ty,
		typename _Ret = conditional_t<_Is_cutoff, index_t, typename _Ty::value_t>>
		MATRICE_GLOBAL_INL auto range(const _Ty& _Pos)->tuple<_Ret, _Ret, _Ret, _Ret> {
		using tuple_type = tuple<_Ret, _Ret, _Ret, _Ret>;
		if constexpr (_Is_cutoff) {
			const auto x = floor<_Ret>(_Pos.x), y = floor<_Ret>(_Pos.y);
			return tuple_type(x - _Radius, x + _Radius + 1, y - _Radius, y + _Radius + 1);
		}
		return tuple_type(_Pos.x - _Radius, _Pos.x + _Radius + 1, _Pos.y - _Radius, _Pos.y + _Radius + 1);
	}

	/**
	 * \brief Check if point '_Pos' in the range '_Shape'.
	 */
	template<typename _Pty, typename _Sty>
	MATRICE_HOST_INL bool range_check(const _Pty& _Pos, const _Sty& _Shape) {
		auto _TL = _Pos - _Radius; 
		auto _RB = _Pos + _Radius;
		if (floor(_TL(0)) < 0 || floor(_TL(1)) < 0 ||
			_RB(0) >= get<0>(_Shape) || _RB(1) >= get<1>(_Shape))
			return false;
		return true;
	}
};

/// <summary>
/// \brief CLASS TEMPLATE, correlation optimizer base
/// </summary>
/// <typeparam name="_Derived"></typeparam>
template<typename _Derived> class _Corr_optim_base {
	using _Myt = _Corr_optim_base;
	using _Mydt = _Derived;
protected:
	using _Mytraits = _Corr_solver_traits<_Mydt>;
public:

	// \brief Number of parameters to be estimated.
	static constexpr auto npar = conditional_size_v<
		_Mytraits::order == 0, 2, _Mytraits::order*6>;
	using value_type = typename _Mytraits::value_type;
	using matrix_type = Matrix<value_type>;
	using matrix_fixed = Matrix_<value_type, npar, npar>;
	using point_type = Vec2_<value_type>;
	using param_type = Vec_<value_type, npar>;
	using jacob_type = Matrix_<value_type, ::dynamic, npar>;
	using vector_type = Matrix_<value_type, ::dynamic, 1>;
	using options_type = _Correlation_options;
	using interp_type = typename _Mytraits::interpolator;
	// \brief Smooth image type supported by spline interpolation.
	using smooth_image_t = interp_type;
	using update_strategy = typename _Mytraits::update_strategy;
	using rect_type = rect<size_t>;

	// \brief Loss function definition.
	struct loss_fn {
		value_type rho(value_type x) noexcept {
			//return scale * (1 - exp(-x/ scale));
			return scale * (1 - 1/(1+x/scale));
			//return x;
		}
		value_type phi(value_type x) noexcept {
			//return exp(-x/scale);
			return 1/sq(1+x/scale);
			//return 1;
		}
		value_type scale = sq(0.001);
	};

	/**
	 *\brief This module is image-wise thread-safe. It allows us to eval. relative deformation of each point between current and reference state.
	 *\param _Ref reference interpolated image;
	 *\param _Cur current interpolated image;
	 *\param _Opt options for the optimizer.
	 */
	_Corr_optim_base(
		const smooth_image_t& _Ref,
		const smooth_image_t& _Cur,
		const options_type& _Opt
	) : _Myopt(_Opt),
		_Myimref(new interp_type(_Ref)),
		_Myimcur(new interp_type(_Cur)),
		_Mysize(_Opt._Radius<<1|1),
		_Myjaco(sq(_Mysize)), _Mydiff(sq(_Mysize)),
		// for robust estimation, Dec/30/2020
		_Myweight(sq(_Mysize)), _Myjaco_tw(sq(_Mysize)) {
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
		DGELOM_CHECK(!_Myref.empty, 
			"Call init(...) before calling method guess(...).");
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

	// \brief for robust estimation, Dec/30/2020
	MATRICE_HOST_INL void set_loss_scale(value_type s) noexcept {
		_Myloss.scale = sq(s);
	}
	// \brief for robust estimation, Dec/30/2020
	// \return [squared_loss, rho_loss, norm2_of_pars] = robust_sol(pars);
	MATRICE_HOST_INL auto robust_sol(param_type& pars);

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
	options_type _Myopt;
	point_type   _Mypos;
	matrix_type  _Myref, _Mycur;
	vector_type  _Mydiff;
	jacob_type   _Myjaco;
	matrix_fixed _Myhess;

	// for robust estimation, Dec/30/2020
	vector_type  _Myweight;
	loss_fn      _Myloss;
	value_type   _Myissd;
	Matrix_<value_type, npar, ::dynamic> _Myjaco_tw;

	// reconstructed reference image with a specified interpolator
	shared_ptr<smooth_image_t> _Myimref;
	// reconstructed current image with a specified interpolator
	shared_ptr<smooth_image_t> _Myimcur;
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
	using typename _Mybase::jacob_type;

	_Corr_solver_impl(const _Myt& _Other) = delete;
	_Corr_solver_impl(_Myt&& _Other) = delete;
	template<typename... _Args>
	_Corr_solver_impl(const _Args&..._args) 
		: _Mybase(_args...) {}

	/**
	 *\brief Evaluate Jacobian contribution at ref. patch.
	 */
	MATRICE_HOST_INL auto& _Diff();

	/**
	 *\brief Construct current image patch.
	 *\param [_Par] displacement {u, v}.
	 */
	MATRICE_HOST_INL auto& _Warp(const param_type& _Par);

protected:
	using typename _Mybase::update_strategy;

	using _Mybase::_Myimref;
	using _Mybase::_Myimcur;
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
	using typename _Mybase::jacob_type;
	using typename _Mybase::matrix_fixed;

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
	 *\brief Thread-safe version of method _Diff().
	 *\param "x, y" The coordinates of a computation point.
	 *\returns Jacobian matrix.
	 */
	MATRICE_HOST_INL auto _Diff(value_type x, value_type y);

	/**
	 *\brief Construct current image patch.
	 *\param [_Par] {u, dudx, dudy, v, dvdx, dvdy}
	 */
	MATRICE_HOST_INL auto& _Warp(const param_type& _Par);

protected:
	using typename _Mybase::update_strategy;

	using _Mybase::_Myimref;
	using _Mybase::_Myimcur;
	using _Mybase::_Mycur;
	using _Mybase::_Myopt;
	using _Mybase::_Mypos;
};

_DETAIL_END
MATRICE_ALG_END(corr)

#include "inline/_optim.inl"