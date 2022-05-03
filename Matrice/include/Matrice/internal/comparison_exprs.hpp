/**********************************************************************
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
**********************************************************************/
#pragma once

#include "expr_base.hpp"
#include "forward.hpp"

MATRICE_NAMESPACE_BEGIN(xpr)

//\brief Less than '<'
struct _Op_lt : public __op__ {
	using retval_t = uint8_t;
	template<typename _Ty>
	MATRICE_GLOBAL_INL retval_t operator()(const _Ty& l, const _Ty& r) const noexcept {
		return l < r;
	}
	template<typename _Ty, typename _Uy>
	MATRICE_GLOBAL_INL retval_t operator()(const _Ty& l, const _Uy& r) const noexcept {
		return l < r;
	}
};
//\brief Less than or equal to '<='
struct _Op_le : public __op__ {
	using retval_t = uint8_t;
	template<typename _Ty, typename _Uy=_Ty>
	MATRICE_GLOBAL_INL retval_t operator()(const _Ty& l, const _Uy& r) const noexcept {
		return l <= r;
	}
};
//\brief Greater than '>'
struct _Op_gt : public __op__ {
	using retval_t = uint8_t;
	template<typename _Ty, typename _Uy=_Ty>
	MATRICE_GLOBAL_INL retval_t operator()(const _Ty& l, const _Uy& r) const noexcept {
		return l > r;
	}
};
//\brief Greater than or equal to '>='
struct _Op_ge : public __op__ {
	using retval_t = uint8_t;
	template<typename _Ty, typename _Uy=_Ty>
	MATRICE_GLOBAL_INL retval_t operator()(const _Ty& l, const _Uy& r) const noexcept {
		return l >= r;
	}
};
//\brief Equal to '=='
struct _Op_eq : public __op__ {
	using retval_t = uint8_t;
	template<typename _Ty, typename _Uy=_Ty>
	MATRICE_GLOBAL_INL retval_t operator()(const _Ty& l, const _Uy& r) const noexcept {
		return l == r;
	}
};
//\brief Not equal to '!='
struct _Op_neq : public __op__ {
	using retval_t = uint8_t;
	template<typename _Ty, typename _Uy=_Ty>
	MATRICE_GLOBAL_INL retval_t operator()(const _Ty& l, const _Uy& r) const noexcept {
		return l != r;
	}
};

template<class _Lhs, class _Rhs, Operator _Op>
class EwiseComparison : public Base<EwiseComparison<_Lhs, _Rhs, _Op>>{
	using _Myt = EwiseComparison;
	using _Mybase = Base<EwiseComparison<_Lhs, _Rhs, _Op>>;
	using _Mylhs_t = conditional_t<is_scalar_v<_Lhs>, Scalar<_Lhs>, const _Lhs&>;
	using _Myrhs_t = conditional_t<is_scalar_v<_Rhs>, Scalar<_Rhs>, const _Rhs&>;
public:
	MATRICE_GLOBAL_FINL EwiseComparison(const _Lhs& lhs, const _Rhs& rhs) noexcept
		:_Mylhs(lhs), _Myrhs(rhs) {
	}

	MATRICE_GLOBAL_FINL auto operator()(size_t idx) const noexcept {
		return _Myop(_Mylhs(idx), _Myrhs(idx));
	}

	MATRICE_GLOBAL_FINL decltype(auto) lhs() const noexcept {
		return (_Mylhs);
	}
	MATRICE_GLOBAL_FINL decltype(auto) rhs() const noexcept {
		return (_Myrhs);
	}
	MATRICE_GLOBAL_FINL auto shape() const noexcept {
		if constexpr (is_scalar_v<_Lhs> || is_scalar_v<_Rhs>) {
			return max(_Myrhs.shape(), _Mylhs.shape());
		}
		else {
			return min(_Myrhs.shape(), _Mylhs.shape());
		}
	}

private:
	_Op _Myop;
	_Mylhs_t _Mylhs;
	_Myrhs_t _Myrhs;
};

template<class T, class U, Operator _Op>
struct traits<EwiseComparison<T, U, _Op>> {
	using value_type = common_type_t<
		typename conditional_t<is_scalar_v<T>, Scalar<T>, T>::value_t, 
		typename conditional_t<is_scalar_v<U>, Scalar<U>, U>::value_t>;
	using retval_type = typename _Op::retval_t;
	static constexpr auto rows = xpr::common_rows_v<T, U>;
	static constexpr auto cols = xpr::common_cols_v<T, U>;
	using return_type = conditional_t<is_scalar_v<T>&&is_scalar_v<U>, 
		Scalar<retval_type>, Matrix_<retval_type, rows, cols>>;
};

MATRICE_NAMESPACE_END(xpr)

DGE_MATRICE_BEGIN
template<typename _Lhs, typename _Rhs, typename _Op = xpr::_Op_lt>
MATRICE_GLOBAL_FINL auto less(const _Lhs& _left, const _Rhs& _right) {
	return xpr::EwiseComparison<_Lhs, _Rhs, _Op>(_left, _right);
}
template<typename _Lhs, typename _Rhs, typename _Op = xpr::_Op_le>
MATRICE_GLOBAL_FINL auto less_equal(const _Lhs& _left, const _Rhs& _right) {
	return xpr::EwiseComparison<_Lhs, _Rhs, _Op>(_left, _right);
}
template<typename _Lhs, typename _Rhs, typename _Op = xpr::_Op_gt>
MATRICE_GLOBAL_FINL auto greater(const _Lhs& _left, const _Rhs& _right) {
	return xpr::EwiseComparison<_Lhs, _Rhs, _Op>(_left, _right);
}
template<typename _Lhs, typename _Rhs, typename _Op = xpr::_Op_ge>
MATRICE_GLOBAL_FINL auto greater_equal(const _Lhs& _left, const _Rhs& _right) {
	return xpr::EwiseComparison<_Lhs, _Rhs, _Op>(_left, _right);
}
template<typename _Lhs, typename _Rhs, typename _Op = xpr::_Op_eq>
MATRICE_GLOBAL_FINL auto equal(const _Lhs& _left, const _Rhs& _right) {
	return xpr::EwiseComparison<_Lhs, _Rhs, _Op>(_left, _right);
}
template<typename _Lhs, typename _Rhs, typename _Op = xpr::_Op_neq>
MATRICE_GLOBAL_FINL auto not_equal(const _Lhs& _left, const _Rhs& _right) {
	return xpr::EwiseComparison<_Lhs, _Rhs, _Op>(_left, _right);
}
DGE_MATRICE_END
