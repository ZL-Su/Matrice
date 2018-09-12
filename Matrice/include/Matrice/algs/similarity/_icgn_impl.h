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
#include "_similarity_traits.h"

MATRICE_ALGS_BEGIN
template<size_t _Options> struct _Options_similarity
{
	// \encapsulate interpolation type, similarity-metric type, etc.
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
namespace details {
template<typename _Derived> class _GaussNewton_base
{
	using derived_type = _Derived;
	using value_type = typename internals::gn_solver_traits<derived_type>::type;
protected:
	using options_type = _Options_similarity<internals::gn_solver_traits<derived_type>::options>;
	enum {dof = internals::static_size<options_type::options>::value};
	using stack_vector = types::Vec_<value_type, dof>;
	using stack_matrix = Matrix_<value_type, dof, dof>;
	using matrix_type  = Matrix<value_type>;
	using const_matrix_reference = const std::add_lvalue_reference_t<matrix_type>;
	struct status_type
	{
		bool _Is_success = false;
		value_type _Value_1 = std::numeric_limits<value_type>::quiet_NaN();
		value_type _Value_2 = std::numeric_limits<value_type>::quiet_NaN();
	};
public:
	using value_t = typename matrix_type::value_t;
	using pointer = std::add_pointer_t<value_t>;
	using point_t = types::Vec_<value_t, internals::static_size<2>::value>;

	_GaussNewton_base(
		const_matrix_reference _Ref, const_matrix_reference _Q,
		const point_t& _Initpos, options_type& _Opts)
		:m_reference(_Ref), m_coeff(_Q), 
		 m_pos(_Initpos), m_options(_Opts){
		m_ksize = m_options() << 1 | 1;
		m_current(m_ksize, m_ksize, zero<value_type>::value);
	}

protected:
	// \engine for iterative solver
	template<typename... _Args> MATRICE_HOST_FINL auto _Solve(_Args... _Args) {
		return static_cast<derived_type*>(this)->_Solver_impl(_Args...);
	}

	// \feature position: m_pos[new] <-- m_pos[old] + delta_pos
	point_t m_pos;

	// \options for iterative refinement
	options_type m_options;

	// \kernel size
	std::size_t m_ksize;

	// \current image buffer: G
	matrix_type m_current;

	// \view of reference image: F
	const_matrix_reference m_reference;

	// \precomputed interpolation coeff.
	const_matrix_reference m_coeff;
};

template<typename _Ty, size_t _Options>
class _GaussNewton_impl_ic MATRICE_NONHERITABLE
	: public _GaussNewton_base<_GaussNewton_impl_ic<_Ty, _Options>>
{
	using base_type = _GaussNewton_base<_GaussNewton_impl_ic<_Ty, _Options>>;
	using value_type = typename base_type::value_t;
	using param_type = types::Vec_(value_type, internals::static_size<6>::value);
public:
	enum { options = _Options };
	using value_t = value_type;
	template<typename... _Args> MATRICE_HOST_FINL
	_GaussNewton_impl_ic(const _Args&... _Args): base_type(_Args...) {}

	MATRICE_HOST_FINL auto _Solver_impl(param_type& _Params);

private:
	value_type m_favg, m_fssd;
	using base_type::m_coeff;
	using base_type::m_ksize;
};
}
MATRICE_ALGS_END
#include "inline\_icgn_impl.inl"
