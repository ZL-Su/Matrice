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
#include "../_simd_accessors.h"

MATRICE_ARCH_BEGIN
_DETAIL_BEGIN
namespace internal{
// Specialization of simd "load" instructions
template<> MATRICE_HOST_FINL
auto _Load(packed_vector<int64_t, 2>, int64_t _Val1, int64_t _Val2) {
	return _mm_set_epi64x(_Val2, _Val1);
}
template<> MATRICE_HOST_FINL
auto _Load(packed_vector<int64_t, 2>, const int64_t* _Src) {
	return _Load(packed_vector<int64_t, 2>{}, _Src[0], _Src[1]);
}
template<> MATRICE_HOST_FINL
auto _Load(packed_vector<double, 2>, const double* _Src) {
	return _mm_load_pd(_Src);
}
template<> MATRICE_HOST_FINL
auto _Load(packed_vector<double, 2>, double _Val1, double _Val2) {
	return _mm_set_pd(_Val2, _Val1);
}
template<> MATRICE_HOST_FINL
auto _Load(packed_vector<float, 4>, float _Src) {
	return _mm_set1_ps(_Src);
}
template<> MATRICE_HOST_FINL
auto _Load(packed_vector<float, 4>, const float* _Src) {
	return _mm_load_ps(_Src);
}
template<> MATRICE_HOST_FINL
auto _Load(packed_vector<float, 4>, float _Val1, float _Val2, float _Val3, float _Val4) {
	return _mm_set_ps(_Val4, _Val3, _Val2, _Val1);
}

// Specialization of simd "set" instructions
template<> MATRICE_HOST_FINL
auto _Set_zero(packed_vector<int64_t, 2>) {
	return _mm_setzero_si128();
}
template<> MATRICE_HOST_FINL
auto _Set_zero(packed_vector<double, 2>) {
	return _mm_setzero_pd();
}
template<> MATRICE_HOST_FINL
auto _Set_zero(packed_vector<float, 4>) {
	return _mm_setzero_ps();
}
}

_DETAIL_END
MATRICE_ARCH_END

#endif