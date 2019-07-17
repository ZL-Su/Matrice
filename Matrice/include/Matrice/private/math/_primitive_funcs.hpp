/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include <cmath>
#if (defined MATRICE_ENABLE_CUDA && !defined __disable_cuda__)
#include <device_functions.h>
#include <thrust\complex.h>
#endif
#include "_config.h"

DGE_MATRICE_BEGIN
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_scalar_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret add(const T1& a, const T2& b) { 
	return a + b; 
}
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_scalar_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret sub(const T1& a, const T2& b) { 
	return a - b; 
}
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_scalar_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret mul(const T1& a, const T2& b) { 
	return a * b; 
}
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_scalar_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret div(const T1& a, const T2& b) { 
	return a / b; 
}
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_scalar_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret max(const T1& a, const T2& b) { 
	return a < b ? b : a; 
}
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_scalar_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret min(const T1& a, const T2& b) { 
	return a < b ? a : b; 
}
template<typename T, typename _Ret = T, typename = enable_if_t<is_scalar_v<T>>>
MATRICE_GLOBAL_FINL constexpr _Ret sqrt(const T& x) { 
	return std::sqrt(x); 
}
template<typename T, typename _Ret = T, typename = enable_if_t<is_scalar_v<T>>>
MATRICE_HOST_FINL constexpr _Ret abs(const T& x) { 
	return std::abs(x); 
}
template<typename T, typename _Ret = T, typename = enable_if_t<is_scalar_v<T>>>
MATRICE_HOST_FINL constexpr _Ret exp(const T& x) { 
	return std::exp(x); 
}
template<typename T, typename _Ret = T, typename = enable_if_t<is_scalar_v<T>>>
MATRICE_HOST_FINL constexpr _Ret log(const T& x) { 
	return std::log(x); 
}
template<typename T, typename _Ret = T, typename = enable_if_t<is_scalar_v<T>>>
MATRICE_HOST_FINL constexpr _Ret log2(const T& x) { 
	return std::log2(x); 
}
template<typename T, typename _Ret = T, typename = enable_if_t<is_scalar_v<T>>>
MATRICE_HOST_FINL constexpr _Ret log10(const T& x) { 
	return std::log10(x); 
}
template<typename T, typename _Ret = T, typename = enable_if_t<is_scalar_v<T>>>
MATRICE_HOST_FINL constexpr _Ret floor(const T& x) { 
	return static_cast<_Ret>(std::floor(x)); 
}
template<typename T, typename _Ret = T, typename = enable_if_t<is_scalar_v<T>>>
MATRICE_HOST_FINL constexpr _Ret ceil(const T& x) { 
	return static_cast<_Ret>(std::ceil(x)); 
}
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_scalar_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret pow(const T1& x, const T2& y) { 
	return std::pow(x,y); 
}
template<typename T1, typename T2, typename _Ret = common_type_t<T1, T2>, typename = enable_if_t<is_scalar_v<_Ret>>>
MATRICE_GLOBAL_FINL constexpr _Ret safe_div(const T1& a, const T2& b) { 
	return b == T2(0) ? _Ret(0) : div(a, b);
}

///<brief> trigonometric functions </brief>
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T sin(const T& x) noexcept {
	return std::sin(x);
}
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T cos(const T& x) noexcept {
	return std::cos(x);
}
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T tan(const T& x) noexcept {
	return std::tan(x);
}
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T arcsin(const T& x) noexcept {
	return std::asin(x);
}
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T arccos(const T& x) noexcept {
	return std::acos(x);
}
template<typename T, MATRICE_ENABLE_IF(is_scalar_v<T>)>
MATRICE_GLOBAL_FINL constexpr T arctan(const T& x) noexcept {
	return std::atan(x);
}

DGE_MATRICE_END