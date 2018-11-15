#pragma once
/**************************************************************************
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
**************************************************************************/
#pragma once
#include "../_matrix_exp.hpp"
#include "../_range.h"

DGE_MATRICE_BEGIN namespace detail {
struct _Tensor_exp_op {
	template<typename T, typename U>
	using _Ewise_mmul = Expr::MatBinaryExpr<T, U, _Exp_op::_Mat_mul<common_value_t<T, U>>>;
};

template<typename _Derived> class _Tensor_exp_base {
	using _Mytraits = tensor_traits<_Derived>;
public:
	MATRICE_GLOBAL_INL auto reduce() const {
		auto _Derived_this = static_cast<const _Derived*>(this);
		auto _Ret = (*_Derived_this)(0).eval();
		for (const auto _Idx : range(1, M*K))
			_Ret = _Ret + (*_Derived_this)(_Idx);
		return (_Ret);
	}
protected:
	std::size_t M, K, N;
};

template<typename T, typename U, typename _Op,
	typename = std::enable_if_t<is_tensor_v<T>&&is_tensor_v<U>>>
class _Tensor_exp : public _Tensor_exp_base<_Tensor_exp<T,U,_Op>> {
	using _Mybase = _Tensor_exp_base<_Tensor_exp<T, U, _Op>>;
public:
	_Tensor_exp(const T& _Left, const U& _Right) 
		: _LHS(_Left), _RHS(_Right) {
		_Mybase::M = _LHS.rows();
		_Mybase::K = _LHS.cols();
		_Mybase::N = _LHS.cols();
	}

	MATRICE_GLOBAL_INL auto operator() (std::size_t _Idx) const {
		return _Op(_LHS(_Idx), _RHS(_Idx));
	}

private:
	const T& _LHS;
	const U& _RHS;
};

} 
DGE_MATRICE_END