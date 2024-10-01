/**********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
#include <cmath>
#include "_config.h"
#include "private/_type_traits.h"
#ifdef MATRICE_ENABLE_CUDA
#include <device_functions.h>
#include <thrust\complex.h>
#endif

DGE_MATRICE_BEGIN
/// <summary>
/// \brief Sum of a scalar number sequence.
/// </summary>
/// <typeparam name="...Ts"></typeparam>
/// <param name="...args">Must be scalar types.</param>
/// <returns> \e.g. 
/// constexpr auto sum = dgelom::add(1, 2, 3.); // 6.0
/// </returns>
template<typename ...Ts>
MATRICE_GLOBAL_FINL constexpr auto add(const Ts&... args)noexcept {
	return (...+ args);
}

/// <summary>
/// \brief Subtract a scalar sequence from a given value.
/// </summary>
/// <typeparam name="T">Must be a scalar type.</typeparam>
/// <typeparam name="...Ts">Must be scalar types.</typeparam>
/// <param name="a">Minuend.</param>
/// <param name="...args">Subtraend sequence.</param>
/// <returns>\e.g.
/// constexpr auto res = dgelom::sub(2., 1, 2); // -1.0 = 2 - 1 - 2
/// </returns>
template<typename T, typename... Ts, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr auto sub(const T& a, const Ts&... args)noexcept {
	return a - (args + ...);
}

/// <summary>
/// \brief Multiplication of a scalar sequence.
/// </summary>
/// <typeparam name="...Ts">Must be scalar types.</typeparam>
/// <param name="...args">A scalar number sequence.</param>
/// <returns> \e.g.
/// constexpr auto res = dgelom::mul(1, 2, 3.); // 6.0 = 1*2*3
/// </returns>
template<typename... Ts>
MATRICE_GLOBAL_FINL constexpr auto mul(const Ts&... args)noexcept {
	return (args * ...);
}

template<typename T1, typename T2, 
	typename _Ret = common_type_t<T1, T2>, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_GLOBAL_FINL constexpr _Ret div(const T1& a, const T2& b)noexcept {
	return a / b;
}
template<typename T1, typename T2, 
	typename _Ret = common_type_t<T1, T2>, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_GLOBAL_FINL constexpr _Ret max(const T1& a, const T2& b)noexcept {
	return a < b ? b : a;
}
template<typename T1, typename T2, 
	typename _Ret = common_type_t<T1, T2>, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_GLOBAL_FINL constexpr _Ret min(const T1& a, const T2& b)noexcept { 
	return a < b ? a : b;
}
template<typename T, typename _Ret = T, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_GLOBAL_FINL constexpr _Ret sqrt(const T& x)noexcept {
	return MATRICE_STD(sqrt)(T(x));
}
template<typename T, typename _Ret = T, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_HOST_FINL constexpr _Ret abs(const T& x)noexcept {
	return MATRICE_STD(abs)(T(x));
}
template<typename T, typename _Ret = T, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_HOST_FINL constexpr _Ret exp(const T& x)noexcept {
	return MATRICE_STD(exp)(T(x));
}
template<typename T, typename _Ret = T, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_HOST_FINL constexpr _Ret log(const T& x)noexcept {
	return MATRICE_STD(log)(T(x));
}
template<typename T, typename _Ret = T, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_HOST_FINL constexpr _Ret log2(const T& x)noexcept {
	return MATRICE_STD(log2)(T(x));
}
template<typename T, typename _Ret = T, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_HOST_FINL constexpr _Ret log10(const T& x)noexcept {
	return MATRICE_STD(log10)(T(x));
}
template<typename T, typename _Ret = T, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_HOST_FINL constexpr _Ret floor(const T& x)noexcept {
	return static_cast<_Ret>(MATRICE_STD(floor)(T(x)));
}
template<typename T, typename _Ret = T, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_HOST_FINL constexpr _Ret ceil(const T& x)noexcept {
	return static_cast<_Ret>(MATRICE_STD(ceil)(T(x)));
}
template<typename T1, typename T2, 
	typename _Ret = common_type_t<T1, T2>, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_GLOBAL_FINL constexpr _Ret pow(const T1& x, const T2& y)noexcept {
	return MATRICE_STD(pow)(T1(x), T2(y));
}
template<typename T1, typename T2, 
	typename _Ret = common_type_t<T1, T2>, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ret>)>
MATRICE_GLOBAL_FINL constexpr _Ret safe_div(const T1& a, const T2& b, T2 thresh = T2(0))noexcept {
	return fabs(b) <= thresh ? _Ret(0) : div(a, b);
}

///<brief> trigonometric functions </brief>
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T sin(const T& x) noexcept {
	return MATRICE_STD(sin)(T(x));
}
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T cos(const T& x) noexcept {
	return MATRICE_STD(cos)(T(x));
}
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T sinc(const T& x) noexcept {
	return sin(pi<T>*x) / (pi<T>*x);
}

/// <summary>
/// \brief Return sine and cose values of a given scalar number.
/// </summary>
/// <typeparam name="T">Must be a scalar type.</typeparam>
/// <param name="x">Any scalar number.</param>
/// <returns>
/// \e.g. auto [sin_x, cos_x] = dgelom::sin_cos(x);
/// </returns>
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr auto sin_cos(const T& x)noexcept {
	return MATRICE_STD(tuple)(sin(x), cos(x));
}
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T tan(const T& x) noexcept {
	return MATRICE_STD(tan)(T(x));
}
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T arcsin(const T& x) noexcept {
	return MATRICE_STD(asin)(T(x));
}
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T arccos(const T& x) noexcept {
	return MATRICE_STD(acos)(T(x));
}
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T arctan(const T& x) noexcept {
	return MATRICE_STD(atan)(T(x));
}

/// <summary>
/// \brief Sum of squares for an input scalar sequence.
/// (Require C++ 17 support.)
/// </summary>
/// <typeparam name="...Ts">Scalar types</typeparam>
/// <param name="...args">A scalar number sequence, such as 'x, y, ...'.</param>
/// <returns> Sum of square sequence: x^2 + y^2 + ...; 
/// \e.g. auto sum_of_squares = dgelom::sqsum(x, y, ...);
/// </returns>
template<typename ...Ts>
MATRICE_GLOBAL_FINL constexpr auto sqsum(const Ts&... args) noexcept {
	return (... + (args * args));
}

DGE_MATRICE_END