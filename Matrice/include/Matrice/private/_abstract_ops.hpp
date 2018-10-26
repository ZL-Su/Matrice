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
#include <array>
#include "_type_traits.h"
#include "../util/genalgs.h"
#include "../util/utils.h"

DGE_MATRICE_BEGIN

using exprs::Expr;

// element-wise addition
template<
	typename _Lhs, typename _Scalar,
	typename _Op = Expr::EwiseBinaryExpr<_Lhs, _Lhs, _Exp_op::_Ewise_sum<typename _Lhs::value_t>>,
	typename = std::enable_if_t<std::is_scalar_v<_Scalar>>>
MATRICE_GLOBAL_FINL auto operator+ (const _Lhs& _opd, const _Scalar _scalar) { return _Op(_scalar, _opd); }
template<
	typename _Rhs,
	typename _Op = Expr::EwiseBinaryExpr<_Rhs, _Rhs, _Exp_op::_Ewise_sum<typename _Rhs::value_t>>>
MATRICE_GLOBAL_FINL auto operator+ (typename _Rhs::value_t _scalar, const _Rhs& _opd) { return _Op(_scalar, _opd); }
template<
	typename _Lhs, typename _Rhs,
	typename value_t = std::common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_sum<value_t>>>
MATRICE_GLOBAL_FINL auto operator+ (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// element-wise subtraction
template<
	typename _Lhs,
	typename _Op = Expr::EwiseBinaryExpr<_Lhs, _Lhs, _Exp_op::_Ewise_min<typename _Lhs::value_t>>>
MATRICE_GLOBAL_FINL auto operator- (const _Lhs& _opd, typename _Lhs::value_t _scalar) { return _Op(-_scalar, _opd); }
template<
	typename _Rhs,
	typename _Op = Expr::EwiseBinaryExpr<_Rhs, _Rhs, _Exp_op::_Ewise_min<typename _Rhs::value_t>>>
MATRICE_GLOBAL_FINL auto operator- (typename _Rhs::value_t _scalar, const _Rhs& _opd) { return _Op(_scalar, _opd); }
template<
	typename _Lhs, class _Rhs,
	typename value_t = std::common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_min<value_t>>>
MATRICE_GLOBAL_FINL auto operator- (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// element-wise multiplication
template<
	typename _Lhs, typename _Scalar,
	typename _Op = Expr::EwiseBinaryExpr<_Lhs, _Lhs, _Exp_op::_Ewise_mul<typename _Lhs::value_t>>,
	typename = std::enable_if_t<std::is_scalar_v<_Scalar>>>
MATRICE_GLOBAL_FINL auto operator* (const _Lhs& _opd, const _Scalar _scalar) { return _Op(_scalar, _opd); }
template<
	typename _Rhs,
	typename _Op = Expr::EwiseBinaryExpr<_Rhs, _Rhs, _Exp_op::_Ewise_mul<typename _Rhs::value_t>>>
MATRICE_GLOBAL_FINL auto operator* (typename _Rhs::value_t _scalar, const _Rhs& _opd) { return _Op(_scalar, _opd); }
template<
	typename _Lhs, class _Rhs,
	typename value_t = std::common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_mul<value_t>>>
MATRICE_GLOBAL_FINL auto operator* (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// element-wise division
template<
	typename _Lhs, typename _Scalar,
	typename _Op = Expr::EwiseBinaryExpr<_Lhs, _Lhs, _Exp_op::_Ewise_mul<typename _Lhs::value_t>>,
	typename = std::enable_if_t<std::is_scalar_v<_Scalar>>>
MATRICE_GLOBAL_FINL auto operator/ (const _Lhs& _opd, const _Scalar _scalar) { return _Op(1./_scalar, _opd); }
template<
	typename _Rhs,
	typename  _Op = Expr::EwiseBinaryExpr<_Rhs, _Rhs, _Exp_op::_Ewise_div<typename _Rhs::value_t>>>
MATRICE_GLOBAL_FINL auto operator/ (typename _Rhs::value_t _scalar, const _Rhs& _opd) { return _Op(_scalar, _opd); }
template<
	typename _Lhs, typename _Rhs,
	typename value_t = std::common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_div<value_t>>>
MATRICE_GLOBAL_FINL auto operator/ (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }

// element-wise sqrt()
template<
	typename _Rhs,
	typename value_t = typename std::enable_if_t<std::is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_sqrt<value_t>>>
MATRICE_GLOBAL_FINL auto sqrt(const _Rhs& _right) { return _Op(_right); }

// element-wise exp()
template<
	typename _Rhs,
	typename value_t = typename std::enable_if_t<std::is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_exp<value_t>>>
MATRICE_GLOBAL_FINL auto exp(const _Rhs& _right) { return _Op(_right); }

// element-wise log()
template<
	typename _Rhs,
	typename value_t = typename std::enable_if_t<std::is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_log<value_t>>>
MATRICE_GLOBAL_FINL auto log(const _Rhs& _right) { return _Op(_right); }

// element-wise log()
template<
	typename _Rhs,
	typename value_t = typename std::enable_if_t<std::is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_abs<value_t>>>
MATRICE_GLOBAL_FINL auto abs(const _Rhs& _right) { return _Op(_right); }

// determinent expression of square matrix
template<
	typename _Rhs,
	typename value_t = typename std::enable_if_t<std::is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>>
MATRICE_HOST_ONLY value_t det_impl(const _Rhs& a);

// transpose expression
template<
	typename _Rhs,
	typename value_t = typename std::enable_if_t<std::is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>, 
	typename     _Op = Expr::MatUnaryExpr<_Rhs, _Exp_op::_Mat_trp<value_t>>>
MATRICE_GLOBAL_FINL auto transpose(const _Rhs& _right) { return _Op(_right); }

// outer product expression : xy^T
template<
	typename _Lhs, typename _Rhs,
	typename value_t = std::common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>, 
	typename     _Op = Expr::MatBinaryExpr<_Lhs, _Rhs, _Exp_op::_Mat_sprmul<value_t>>>
MATRICE_GLOBAL_FINL auto outer_product(const _Lhs& _left, const _Rhs& _right) {
	return _Op(_left, _right); 
}

_INTERNAL_BEGIN
// helper operators
template<typename _InIt, typename = std::enable_if_t<std::is_pointer_v<_InIt>>>
MATRICE_GLOBAL_FINL void _Conformity_check(_InIt _Left, _InIt _Right) {
#ifdef _DEBUG
	if (_Left != _Right) throw std::runtime_error("Inconsistent iterators!");
#endif
	_Left = _Right;
}
template<typename _InIt, typename = std::enable_if_t<std::is_pointer_v<_InIt>>>
MATRICE_GLOBAL_FINL _InIt _Proxy_checked(const _InIt _Right) {
	return (_Right);
}
template<typename _Valty, size_t _N, typename = std::enable_if_t<std::is_arithmetic_v<_Valty>>>
MATRICE_GLOBAL_FINL auto _Fill_array(const _Valty* _First) {
#ifdef _DEBUG
	if (_First + _N - 1 == nullptr) throw std::runtime_error("Unconformable shape!");
#endif
	std::array<_Valty, _N> _Ret;
	std::_Copy_unchecked(_First, _First + _N, _Ret.data());
	return (_Ret);
}

template<std::size_t _P> struct _Matrix_norm_impl {
	template<typename _Mty>
	MATRICE_GLOBAL_FINL static constexpr auto value(const _Mty& _A) {
		static_assert(is_matrix_v<_Mty>, "_Mty is not matrix type.");

		auto _Ret = typename _Mty::value_type(0);
		for(const auto& _Val : _A){ 
			_Ret += detail::_Powers_n<_P>::value(abs(_Val)); 
		};
		return (pow(_Ret, decltype(_Ret)(1)/_P));
	}
};
template<> struct _Matrix_norm_impl<0> { //infinity norm
	template<typename _Mty>
	MATRICE_GLOBAL_FINL static constexpr auto value(const _Mty& _A) {
		static_assert(is_matrix_v<_Mty>, "_Mty is not matrix type.");

		auto _Ret = typename _Mty::value_type(0);
		for (std::size_t _Idx = 0; _Idx < _A.rows(); ++_Idx) {
			_Ret = max(_Ret, abs(_A.rview(_Idx)).sum());
		}
		return (_Ret);
	}
};
template<> struct _Matrix_norm_impl<1> {
	template<typename _Mty>
	MATRICE_GLOBAL_FINL static constexpr auto value(const _Mty& _A) {
		static_assert(is_matrix_v<_Mty>, "_Mty is not matrix type.");

		auto _Ret = typename _Mty::value_type(0);
		for (std::size_t _Idx = 0; _Idx < _A.cols(); ++_Idx) {
			_Ret = max(_Ret, abs(_A.cview(_Idx)).sum());
		}
		return (_Ret);
	}
};
_INTERNAL_END

DGE_MATRICE_END