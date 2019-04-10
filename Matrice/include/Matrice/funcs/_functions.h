/**************************************************************************
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
**************************************************************************/
#pragma once
#include <type_traits>
#include <exception>
#include "../util/_macros.h"
#include "../util/_std_wrapper.h"
#include "../util/utils.h"
#include "../arch/ixpacket.h"
#include "../private/_type_traits.h"
#include "../private/_size_traits.h"
#include "../private/math/_config.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty> MATRICE_HOST_INL constexpr
_Ty _Det(const _Ty* data, int n, enable_if_t<is_float32_v<_Ty>>* = 0) {
	return fblas::_sdetm(static_cast<float*>(data), n);
}
template<typename _Ty> MATRICE_HOST_INL constexpr
_Ty _Det(const _Ty* data, int n, enable_if_t<is_float64_v<_Ty>>* = 0) {
	return fblas::_ddetm(static_cast<double*>(data), n);
}
_DETAIL_END

/**
 * \SIMD supported determinant of a matrix or its derivatives.
 */
template<typename _T, MATRICE_ENABLE_IF(is_matrix_convertible_v<_T>)>
MATRICE_HOST_INL constexpr typename _T::value_t det(const _T& _x) {
	return detail::_Det(_x.eval().data(), _x.rows());
}

/**
 * \SIMD supported dot-product of two vectors.
 */
template<typename _T, typename _U, 
	typename _Vty = common_type_t<typename _T::value_t, typename _U::value_t>>
MATRICE_HOST_INL constexpr _Vty dot(const _T& _x, const _U& _y) {
#ifdef _DEBUG
	DGELOM_CHECK(_x.size()== _y.size(), "Oops, non-consistent size error.");
#endif

	auto _Ret = zero<_Vty>;

#if   MATRICE_MATH_KERNEL == MATRICE_USE_MKL

#elif MATRICE_MATH_KERNEL == MATRICE_USE_FKL

#else //MATRICE_MATH_KERNEL == MATRICE_USE_NAT
	using packet_type = simd::Packet_<_Vty>;

	const auto _Left = _x.eval();
	const auto _Right = _y.eval();
	const auto _Size = min(_Left.size(), _Right.size());
	const auto _Vsize = simd::vsize<packet_type::size>(_Size);
	for (auto _Idx = 0; _Idx < _Vsize; _Idx += packet_type::size) {
		_Ret += (
			packet_type(_Left.data() + _Idx)*
			packet_type(_Right.data() + _Idx)).reduce();
	}
	for (auto _Idx = _Vsize; _Idx < _Size; ++_Idx) {
		_Ret += _Left(_Idx)*_Right(_Idx);
	}
#endif

	return (_Ret);
}

/**
 * 1D gaussian kernel function: 
 *		g(x;\mu,\sigma) = \frac{1}{\sigma \sqrt{2\pi}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}
 */
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
MATRICE_HOST_INL constexpr auto gaussian(_Ty x, _Ty _M, _Ty _S) {
	return (0.3989422804*exp(-0.5*pow((x - _M) / _S, 2)) / _S);
}

/**
 * 2D gaussian kernel function with uniform Mean and STD for both directions
 */
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
MATRICE_HOST_INL constexpr auto gaussian(_Ty x, _Ty y, _Ty _M, _Ty _S) {
	const auto _S2 = pow(_S, 2);
	return (0.1591549431*exp(-0.5*(pow(x - _M, 2) + pow(y - _M, 2)) / _S2)/ _S2);
}
DGE_MATRICE_END
