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
#pragma once
#include "../../core/matrix.h"
#include "../../core/vector.h"
#include "../../core/solver.h"
#include "../../core/tensor.h"
#include "../../private/math/_linear.h"
#include "../imageproc.hpp"
#include "_similarity_traits.h"

MATRICE_ALGS_BEGIN
template<size_t _Options> struct _Iterative_conv_options
{
	// \interpolation type
	enum { options = _Options };

	// \subset radius, default value is 10
	size_t _My_radius = 10;

	// \spacing between adjacent nodes, valid only for regular lattice
	size_t _My_stride =  7;

	// \maximum iterations, default is 50
	size_t _My_maxits = 50;

	// \termination tolerance condition of refinement iteration
	default_type _My_epsilon = 1.0E-6;

	MATRICE_GLOBAL_FINL auto& operator()() { return (_My_radius); }
	MATRICE_GLOBAL_FINL const auto& operator()() const { return (_My_radius); }
	MATRICE_GLOBAL_FINL operator default_type() { return _My_epsilon; }
};

namespace detail {
template<typename _Itptag> struct _Conv_border_size {};
template<> struct _Conv_border_size<_TAG bicspl_tag> {
	static constexpr auto lower = 1, upper = 2;
};
template<> struct _Conv_border_size<_TAG biqspl_tag> {
	static constexpr auto lower = 2, upper = 3;
};
template<> struct _Conv_border_size<_TAG bisspl_tag> {
	static constexpr auto lower = 3, upper = 4;
};

// \TEMPLATE base class for Gaussian-Newton algorithm [thread-safe]
template<typename _Derived> class _Iterative_conv_base
{
	using _Mydt = _Derived;
	using _Mytraits = internal::conv_solver_traits<_Mydt>;
protected:
	enum {DOF = _Mytraits::order*6}; //DOF
	using value_type = typename _Mytraits::value_type;
	using stack_vector = Vec_<value_type, DOF>;
	using stack_matrix = Matrix_<value_type, DOF, DOF>;
	using matrix_type  = Matrix<value_type>;
	using const_matrix_reference = const std::add_lvalue_reference_t<matrix_type>;
	using param_type = stack_vector;
	using linear_solver_type = matrix_decomp<stack_matrix, _TAG _Linear_spd_tag>;
	using options_type = _Iterative_conv_options<_Mytraits::interp>;
	using interp_category = typename _Mytraits::interp_category;
	using interp_type = typename interpolation<value_type, interp_category>::type;

	struct status_type {
		bool _Is_success = false;
		value_type _Value_1 = std::numeric_limits<value_type>::quiet_NaN();
		value_type _Value_2 = std::numeric_limits<value_type>::quiet_NaN();
	};
public:
	using value_t = typename matrix_type::value_t;
	using pointer = std::add_pointer_t<value_t>;
	using point_t = Vec_<value_t, internal::static_size<2>::value>;

	_Iterative_conv_base(
		const multi_matrix<value_type>& _F,
		const interp_type& _Itp,
		const point_t& _Pos,
		const options_type& _Opts)
		:m_reference(_F), _Myitp(std::make_shared<interp_type>(_Itp)),
		m_pos(_Pos), m_options(_Opts),
		m_ksize(_Opts() << 1 | 1), _Mysolver(_Myhess) {
		_Myhess.format = symm;
		_Myref.create(m_ksize, m_ksize);
		_Mycur.create(m_ksize, m_ksize);
	}
	_Iterative_conv_base(
		const multi_matrix<value_type>& _F,
		const std::shared_ptr<interp_type>& _Itp,
		const point_t& _Pos,
		const options_type& _Opts)
		:m_reference(_F), _Myitp(_Itp), 
		m_pos(_Pos), m_options(_Opts),
		m_ksize(_Opts() << 1 | 1), _Mysolver(_Myhess) {
		_Myhess.format = symm;
		_Myref.create(m_ksize, m_ksize);
		_Mycur.create(m_ksize, m_ksize);
	}

	// \set initial pos
	MATRICE_HOST_FINL auto& pos() { return (m_pos); }
	// \get refined pos
	MATRICE_HOST_FINL const auto& pos() const { return (m_pos); }

protected:
	// \engine for iterative solver
	template<typename... _Args> 
	MATRICE_HOST_FINL auto _Solve(_Args... _Args) {
		return static_cast<_Mydt*>(this)->_Impl(_Args...);
	}

	MATRICE_HOST_FINL auto _Update_subset(const param_type& _P);

	// \feature position: m_pos[new] <-- m_pos[old] + delta_pos
	point_t m_pos;

	// \options for iterative refinement
	options_type m_options;

	// \kernel size
	std::size_t m_ksize;

	// \reference and current image patch buffers: f and g
	matrix_type _Myref, _Mycur;

	// \view of reference image and its gradients: F, dFdx, dFdy
	const multi_matrix<value_type>& m_reference;

	// \precomputed interpolation coeff.
	std::shared_ptr<interp_type> _Myitp;

	// \Jacobian
	tensor_<value_type, 1, DOF> _Myjaco;
	// \Hessian
	stack_matrix _Myhess;

	linear_solver_type _Mysolver;
};

// TEMPLATE impl class for IC-GN optimization [thread-safe]
template<typename _Ty = float, 
	typename _Tag = _TAG bicspl_tag, std::size_t _Ord = 1>
class _Invcomp_conv_impl MATRICE_NONHERITABLE
	: public _Iterative_conv_base<_Invcomp_conv_impl<_Ty, _Tag, _Ord>>
{
	using _Myt = _Invcomp_conv_impl;
	using _Mybase = _Iterative_conv_base<_Myt>;
	using typename _Mybase::options_type;
	using typename _Mybase::const_matrix_reference;
	using typename _Mybase::param_type;
	using typename _Mybase::stack_vector;
	using typename _Mybase::interp_type;
	using value_type = typename _Mybase::value_t;
public:
	enum { options = options_type::options };
	using value_t = value_type;
	using options_t = options_type;
	using param_t = param_type;
	using typename _Mybase::point_t;

	/**
	 * _Ref : reference image package, includes reference image and its gradients w.r.t. x and y
	 * _Itp : interpolator of current image
	 * _Pos : the anchor position of the reference patch
	 * _Opts: solver settings
	 */
	MATRICE_HOST_FINL _Invcomp_conv_impl(
		const multi_matrix<value_t>& _Ref, 
		const std::shared_ptr<interp_type>& _Itp, 
		point_t _Pos, options_t _Opts = options_t()) noexcept
		: _Mybase(_Ref, _Itp, _Pos, _Opts) {
		_Init();
	}
	MATRICE_HOST_FINL _Invcomp_conv_impl(
		const multi_matrix<value_t>& _Ref,
		const interp_type& _Itp,
		point_t _Pos, options_t _Opts = options_t()) noexcept
		: _Mybase(_Ref, _Itp, _Pos, _Opts) {
		_Init();
	}

	MATRICE_HOST_FINL auto _Impl(param_type& _Params);

private:
	MATRICE_HOST_FINL auto _Init();
	MATRICE_HOST_FINL auto _Update();
};
}
MATRICE_ALGS_END
#include "inline\_iterative_conv_impl.inl"
