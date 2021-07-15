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
#include "_simd_traits.hpp"

MATRICE_ARCH_BEGIN
_DETAIL_BEGIN

template<typename... _Args>
MATRICE_HOST_FINL auto _Load(_Args... _) {};

template<packetable_scalar _Ty,
	uint8_t _N = packet_traits<_Ty>::size> struct _Accessor {};

template<packetable_scalar _Ty>
struct _Accessor<_Ty, 2> : packed_vector<_Ty, 2>
{
	using _Mybase = packed_vector<_Ty, 2>;

	template<typename... _Args>
	MATRICE_HOST_FINL auto operator()(_Args&&... _Src) noexcept {
		return _Load(_Mybase{}, _Src...);
	}
};

template<> MATRICE_HOST_FINL 
auto _Load(packed_vector<double, 2>, const double* _Src) {
	return _mm_load_pd(_Src);
}
template<> MATRICE_HOST_FINL 
auto _Load(packed_vector<double, 2>, double _Val1, double _Val2) {
	return _mm_set_pd(_Val2, _Val1);
}

_DETAIL_END

template<typename _Ty, uint8_t _Size= packet_traits<_Ty>::size>
using accessor = detail::_Accessor<_Ty, _Size>;

MATRICE_ARCH_END

#endif // MATRICE_SIMD_ARCH
