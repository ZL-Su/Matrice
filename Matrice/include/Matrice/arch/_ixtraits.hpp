/*********************************************************************
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
**********************************************************************/
#pragma once
#include "private/_type_traits.h"
#include "util/_macros.h"

#ifdef MATRICE_SIMD_ARCH
#include <mmintrin.h>    //_m64
#include <emmintrin.h>   //_m128
#include <immintrin.h>   //_m256
#include <zmmintrin.h>   //_m512
MATRICE_ARCH_BEGIN
#pragma region <!-- simd traits -->
template<typename T, int _Elems> struct conditional {
	static_assert(true, "Invalid _Elems in simd::conditional<T, _Elems>.");
};
template<typename T> struct conditional<T, 2> {
	using type = dgelom::conditional_t<is_common_int64_v<T>, __m128i, dgelom::conditional_t<is_float32_v<T>, __m64, dgelom::conditional_t<is_float64_v<T>, __m128d, void>>>;
};
template<typename T> struct conditional<T, 4> {
	using type = dgelom::conditional_t<is_common_int64_v<T>, __m256i, dgelom::conditional_t<is_float32_v<T>, __m128, dgelom::conditional_t<is_float64_v<T>, __m256d, void>>>;
};
template<typename T> struct conditional<T, 8> {
	using type = dgelom::conditional_t<is_common_int64_v<T>, __m512i, dgelom::conditional_t<is_float32_v<T>, __m256, dgelom::conditional_t<is_float64_v<T>, __m512d, void>>>;
};
template<typename T> struct conditional<T, 16> {
	using type = dgelom::conditional_t<dgelom::is_float32_v<T>, __m512, dgelom::conditional_t<dgelom::is_float64_v<T>, __m512d, void>>;
};
template<typename T, int _Elems> 
using conditional_t = typename conditional<T, _Elems>::type;

template<typename T> struct packet_size {
	static_assert(true, "Unsupported data type T in simd::packet_size<T>.");
};
template<> struct packet_size<float> {
	static constexpr int value =
#if MATRICE_SIMD_ARCH == MATRICE_SIMD_AVX
		1 << MATRICE_SIMD_AVX
#elif MATRICE_SIMD_ARCH == MATRICE_SIMD_AVX512
		1 << MATRICE_SIMD_AVX512
#else
		1 << MATRICE_SIMD_SSE
#endif
		;
};
template<> struct packet_size<double> {
	static constexpr int value =
#if MATRICE_SIMD_ARCH == MATRICE_SIMD_AVX
		1 << ~-MATRICE_SIMD_AVX
#elif MATRICE_SIMD_ARCH == MATRICE_SIMD_AVX512
		1 << ~-MATRICE_SIMD_AVX512
#else
		1 << ~-MATRICE_SIMD_SSE
#endif
		;
};

template<typename T>
MATRICE_HOST_INL constexpr auto packet_size_v = packet_size<T>::value;

template<typename T, int _Elems> struct simd_traits
{
	struct is_sfp8 { enum { value = is_float32<T>::value && _Elems == 8 }; };
	struct is_dfp8 { enum { value = is_float64<T>::value && _Elems == 8 }; };
};
template<typename T, int _Elems> MATRICE_HOST_INL
constexpr auto is_simd_sfx8_v = simd_traits<T, _Elems>::is_sfp8::value;
template<typename T, int _Elems> MATRICE_HOST_INL
constexpr auto is_simd_dfx8_v = simd_traits<T, _Elems>::is_dfp8::value;
#pragma endregion
MATRICE_ARCH_END
#endif // __AVX__