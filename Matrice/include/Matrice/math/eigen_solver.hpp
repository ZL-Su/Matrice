/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2022, Zhilong(Dgelom) Su, all rights reserved.

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
*********************************************************************/
#pragma once

#include "internal/expr_base.hpp"
#include "private/_type_traits.h"
#include "core/matrix.h"

MATRICE_NAMESPACE_BEGIN(detail)

/// <summary>
/// \brief Expression template for building Eigen Problem Solver.
/// The solver expression accepts a matrix or matrix expression as input,
/// and thus it cannot be transferred out the scope where the input lives
/// due to the lifetime problem unless you really know what you do.
/// 
/// As an expression, the return behavior of the eigensolver relies on
/// the lvalue type to be used to hold the result. E.g., the eigenvalues
/// are returned if the solver is assigned to a vector with a size equal
/// to the column of the given matrix. Of course, some methods like 
/// 'eval(...)' and 'compute(...)' are also provided for explicit
/// computation with the specified options (such as if the eigenvectors
/// are required or not).
/// </summary>
/// <typeparam name="_Xtype">
/// Abstract template type, can be Matrix_ or its expression.
/// </typeparam>
template<class _Xtype>
class _Eigen_solver : public xpr::__xpr__ {
	using _Mytraits = conditional_t<is_matrix_v<_Xtype>, 
		matrix_traits<_Xtype>, expression_traits<_Xtype>>;
public:
	using value_type = _Mytraits::value_type;

	_Eigen_solver(const _Xtype& a) noexcept
		:_Myobj(a) {
	}

	MATRICE_HOST_INL auto compute(bool eigenvector = false) {
		auto _Mat = _Myobj.eval();
	}

private:
	const _Xtype& _Myobj;

};

MATRICE_NAMESPACE_END(detail)