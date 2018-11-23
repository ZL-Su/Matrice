/*  *************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
*	*************************************************************************/
#pragma once
#include "../private/math/_linear_solver.h"

DGE_MATRICE_BEGIN
_TYPES_BEGIN

struct Solver MATRICE_NONHERITABLE
{
	template<typename _Op = detail::LinearOp::Auto<Matrix_<default_type, __, __>>> 
	class Linear_ : public detail::SolverBase<Linear_<_Op>>
	{
		using _Mybase = detail::SolverBase<Linear_<_Op>>;
		using value_t = typename _Op::value_t;
		using typename _Mybase::Options;
	public:
		template<typename... _Args> MATRICE_GLOBAL_FINL
		constexpr Linear_(const _Args&... args) : m_op(args...) {};
		template<typename... _Args> MATRICE_GLOBAL_FINL
		constexpr auto solve(const _Args&... args) { return _Mybase::_Impl(args...); }

		template<typename _Rhs>
		MATRICE_GLOBAL_INL auto& operator()(_Rhs& _X) { return m_op(_X); }

		_Op m_op;
	};
};
_TYPES_END

/**
 *\Linear algebra system kernel.
 */
using linear_alg_op = detail::LinearOp;

/**
 *\Linear solver, default _Op is auto-solver-kernel.
 */
template<
	typename _Op = linear_alg_op::Auto<types::Matrix_<default_type, __, __>>,
	typename = std::enable_if_t<is_matrix_v<typename _Op::_Mty>>>
using linear_solver = types::Solver::Linear_<_Op>;

DGE_MATRICE_END