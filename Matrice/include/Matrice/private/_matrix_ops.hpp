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
#include "../util/genalgs.h"
#include "../util/utils.h"

DGE_MATRICE_BEGIN

_INTERNAL_BEGIN
// *\helper operators
template<typename _InIt, MATRICE_ENABLE_IF(is_pointer_v<_InIt>)>
MATRICE_GLOBAL_FINL void _Conformity_check(_InIt _Left, _InIt _Right) {
#ifdef _DEBUG
	DGELOM_CHECK(_Left == _Right, "Inconsistent iterators.");
#endif
	_Left = _Right;
}
template<typename _InIt, MATRICE_ENABLE_IF(is_pointer_v<_InIt>)>
MATRICE_GLOBAL_FINL _InIt _Proxy_checked(const _InIt _Right) {
	return (_Right);
}
template<typename _Valty, size_t _N, 
	MATRICE_ENABLE_IF(is_arithmetic_v<_Valty>)>
MATRICE_GLOBAL_FINL auto _Fill_array(const _Valty* _First) {
#ifdef _DEBUG
	DGELOM_CHECK(_First+_N-1, "The shapes of the source and destination are inconsistent.")
#endif
	std::array<_Valty, _N> _Ret;
	std::_Copy_unchecked(_First, _First + _N, _Ret.data());
	return (_Ret);
}
// *\determinent expression of square matrix
template<typename _Rhs,
	typename value_t = enable_if_t<is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>>
MATRICE_HOST_ONLY value_t det_impl(const _Rhs& a);


/**
 *\brief Evaluate the p-norm for a given matrix or vector.
 *\param [_P] 
 */
template<size_t _P> struct _Matrix_norm_impl {
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
/**
 *\brief Evaluate the 0-norm for a given matrix or vector.
 */
template<> struct _Matrix_norm_impl<0> { //infinity norm
	template<typename _Mty>
	MATRICE_GLOBAL_FINL static constexpr auto value(const _Mty& _A) {
		static_assert(is_matrix_v<_Mty>, "_Mty is not matrix type.");

		auto _Ret = typename _Mty::value_type(0);
		for (size_t _Idx = 0; _Idx < _A.rows(); ++_Idx) {
			_Ret = max(_Ret, abs(_A.rview(_Idx)).sum());
		}
		return (_Ret);
	}
};
/**
 *\brief Evaluate the 1-norm for a given matrix or vector.
 */
template<> struct _Matrix_norm_impl<1> {
	template<typename _Mty>
	MATRICE_GLOBAL_FINL static constexpr auto value(const _Mty& _A) {
		static_assert(is_matrix_v<_Mty>, "_Mty is not matrix type.");

		auto _Ret = typename _Mty::value_type(0);
		for (size_t _Idx = 0; _Idx < _A.cols(); ++_Idx) {
			_Ret = max(_Ret, abs(_A.cview(_Idx)).sum());
		}
		return (_Ret);
	}
};

_INTERNAL_END

DGE_MATRICE_END