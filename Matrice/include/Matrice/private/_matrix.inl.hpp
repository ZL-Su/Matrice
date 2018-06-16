/***************************************************************************
This file is part of Matrice, an effcient and elegant C++ library for SC.
Copyright(C) 2018, Zhilong (Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) 
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/
#pragma once

#include "_expr_type_traits.h"

MATRICE_NAMESPACE_BEGIN_TYPES

using exprs::Expr;
// element-wise addition
template<typename _Rhs,
	      typename _Op = Expr::EwiseBinaryExpr<_Rhs, _Rhs, Expr::Op::EwiseSum<typename _Rhs::value_t>>>
MATRICE_GLOBAL_INL auto operator+ (typename _Rhs::value_t _scalar, const _Rhs& _opd) { return _Op(_scalar, _opd); }
template<typename _Lhs, class _Rhs,
	     typename value_t = std::common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>,
	     typename     _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, Expr::Op::EwiseSum<value_t>>>
MATRICE_GLOBAL_INL auto operator+ (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// element-wise subtraction
template<typename _Rhs,
	      typename _Op = Expr::EwiseBinaryExpr<_Rhs, _Rhs, Expr::Op::EwiseMin<typename _Rhs::value_t>>>
MATRICE_GLOBAL_INL auto operator- (typename _Rhs::value_t _scalar, const _Rhs& _opd) { return _Op(_scalar, _opd); }
template<typename _Lhs, class _Rhs,
	     typename value_t = std::common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>,
	     typename     _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, Expr::Op::EwiseMin<value_t>>>
MATRICE_GLOBAL_INL auto operator- (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// element-wise multiplication
template<typename _Rhs,
	      typename _Op = Expr::EwiseBinaryExpr<_Rhs, _Rhs, Expr::Op::EwiseMul<typename _Rhs::value_t>>>
MATRICE_GLOBAL_INL auto operator* (typename _Rhs::value_t _scalar, const _Rhs& _opd) { return _Op(_scalar, _opd); }
template<typename _Lhs, class _Rhs,
		   typename value_t = std::common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>, 
	      typename     _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, Expr::Op::EwiseMul<value_t>>>
MATRICE_GLOBAL_INL auto operator* (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// element-wise division
template<typename _Rhs, 
	      typename  _Op = Expr::EwiseBinaryExpr<_Rhs, _Rhs, Expr::Op::EwiseDiv<typename _Rhs::value_t>>>
MATRICE_GLOBAL_INL auto operator/ (typename _Rhs::value_t _scalar, const _Rhs& _opd) { return _Op(_scalar, _opd); }
template<typename _Lhs, typename _Rhs, 
	     typename value_t = std::common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>, 
	     typename     _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, Expr::Op::EwiseDiv<value_t>>>
MATRICE_GLOBAL_INL auto operator/ (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }

MATRICE_NAMESPACE_END_TYPES

MATRICE_NAMESPACE_BEGIN_
using exprs::Expr;

// determinent expression of square matrix
template<typename _Rhs,
	     typename value_t = typename std::enable_if<std::is_scalar<typename _Rhs::value_t>::value, typename _Rhs::value_t>::type>
MATRICE_HOST_ONLY value_t det_impl(const _Rhs& a);
// transpose expression
template<typename _Rhs,
	     typename value_t = typename std::enable_if<std::is_scalar<typename _Rhs::value_t>::value, typename _Rhs::value_t>::type, 
	     typename     _Op = Expr::MatUnaryExpr<_Rhs, Expr::Op::MatTrp<value_t>>>
MATRICE_GLOBAL_INL auto transpose(const _Rhs& _right) { return _Op(_right); }
// outer product expression : xy^T
template<typename _Lhs, typename _Rhs,
	     typename value_t = std::common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>, 
	     typename     _Op = Expr::MatBinaryExpr<_Lhs, _Rhs, Expr::Op::SpreadMul<value_t>>>
MATRICE_GLOBAL_INL auto outer_product(const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }

template<typename _InIt, typename = std::enable_if_t<std::is_pointer_v<_InIt>>>
MATRICE_GLOBAL_INL void _Conformity_check(_InIt _Left, _InIt _Right)
{
	if (_Left != _Right) std::runtime_error("Inconsistent iterators!");
	_Left = _Right;
}
template<typename _InIt, typename = std::enable_if_t<std::is_pointer_v<_InIt>>>
MATRICE_GLOBAL_INL _InIt _Proxy_checked(const _InIt _Right)
{
	return (_Right);
}

_MATRICE_NAMESPACE_END