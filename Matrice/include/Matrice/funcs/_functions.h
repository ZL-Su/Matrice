/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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
#include "util/_macros.h"
#include "util/_std_wrapper.h"
#include "util/utils.h"
#include "private/_type_traits.h"
#include "private/_size_traits.h"
#include "private/math/_config.h"
#ifdef MATRICE_SIMD_ARCH
#include "arch/ixpacket.h"
#endif

DGE_MATRICE_BEGIN
#ifdef MATRICE_SIMD_ARCH
/**
 * \SIMD supported dot-product of two vectors.
 */
template<typename _T, typename _U, 
	typename value_t = common_type_t<typename _T::value_t, typename _U::value_t>>
MATRICE_HOST_INL constexpr value_t dot(const _T& _x, const _U& _y) {
#ifdef _DEBUG
	DGELOM_CHECK(_x.size()== _y.size(), "Oops, non-consistent size error.");
#endif
	using packet_type = simd::Packet_<value_t>;

	const auto _Left = _x.eval();
	const auto _Right = _y.eval();
	const auto _Size = min(_Left.size(), _Right.size());
	const auto _Vsize = simd::vsize<packet_type::size>(_Size);

	auto _Ret = zero<value_t>;
	for (auto _Idx = 0; _Idx < _Vsize; _Idx += packet_type::size) {
		_Ret += (packet_type(_Left.data() + _Idx)*
			packet_type(_Right.data() + _Idx)).reduce();
	}
	for (auto _Idx = _Vsize; _Idx < _Size; ++_Idx) {
		_Ret += _Left(_Idx)*_Right(_Idx);
	}

	return (_Ret);
}
#endif

/**
 *\brief 1D gaussian kernel function: 
 *		 g(x;\mu,\sigma) = \frac{1}{\sigma \sqrt{2\pi}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}
 */
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
MATRICE_HOST_INL constexpr _Ty gaussian(_Ty x, _Ty _M, _Ty _S) {
	return (0.3989422804*exp(-0.5*pow((x - _M) / _S, 2)) / _S);
}

/**
 *\brief 2D gaussian kernel function with uniform Mean and STD for both directions
 */
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
MATRICE_HOST_INL constexpr _Ty gaussian(_Ty x, _Ty y, _Ty _M, _Ty _S) {
	const auto _S2 = pow(_S, 2);
	return (0.1591549431*exp(-0.5*(pow(x - _M, 2) + pow(y - _M, 2)) / _S2)/ _S2);
}

/**
 *\brief Clamp an integer x to _N (must be 2^n -1) s.t., for _N > 0, it returns _N if x > _N, or x otherwise; for _N = 0, it returns 0 if x < 0, or x otherwise
 */
template<uint64_t _N>
MATRICE_GLOBAL_FINL constexpr int64_t clamp_to(int64_t x) noexcept {
	return (((_N - x) >> 31) | x) & _N;
}
template<>
MATRICE_GLOBAL_FINL constexpr int64_t clamp_to<0>(int64_t x) noexcept {
	return ((-x) >> 31) & x;
}

/**
 *\brief Fast divide an integer x by 255.
 */
template<typename _Ty, MATRICE_ENABLE_IF(is_integral_v<_Ty>)>
MATRICE_GLOBAL_FINL constexpr _Ty fast_div_255(_Ty x) noexcept {
	return (((x)+(((x)+257) >> 8)) >> 8);
}

/**
 *\brief Fast check if an integer x in range of [min, max].
 */
template<typename _Ty, MATRICE_ENABLE_IF(is_integral_v<_Ty>)>
MATRICE_GLOBAL_FINL constexpr bool fast_range_check(_Ty x, _Ty min, _Ty max) noexcept {
	return ((x - min) | (max - x)) >= 0;
}

DGE_MATRICE_END
