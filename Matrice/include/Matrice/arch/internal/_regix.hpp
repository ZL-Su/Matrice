/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
*********************************************************************/
#pragma once
#ifdef MATRICE_SIMD_ARCH
#include "arch/_simd_traits.hpp"

MATRICE_ARCH_BEGIN
_DETAIL_BEGIN
/// <summary>
/// \brief Define an uniform interface for SIMD intrinsic vector.
/// </summary>
/// <typeparam name="_Intr">Supported intrinsic type.</typeparam>
template<typename _Intr> struct _Regix;
//{ static_assert(true, "Unknown intrinsic vector type in _Regix<_Intr>."); };

/// <summary>
/// \brief Specializations.
/// </summary>
template<> struct _Regix<int32_t> {};
template<> struct _Regix<int64_t> {};
template<> struct _Regix<float_t> {};
template<> struct _Regix<double_t> {};

template<> struct _Regix<__m128> {
	using value_t = float_t;
	using pointer = value_t*;

	MATRICE_HOST_FINL decltype(auto) operator<<(const pointer const _Src) noexcept {
		_Myreg = (_mm_load_ps)(_Src);
		return (*this);
	}
	MATRICE_HOST_FINL auto operator>>(pointer _Dst) const noexcept {
		(_mm_store_ps)(_Dst, _Myreg);
		return _Dst;
	}
	MATRICE_HOST_FINL const auto& operator()() const noexcept {
		return (_Myreg);
	}
	MATRICE_HOST_FINL auto& operator()() noexcept {
		return (_Myreg);
	}

	__m128 _Myreg;
};
_DETAIL_END
MATRICE_ARCH_END

#endif