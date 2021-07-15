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
#include "private/_type_traits.h"
#include "private/_size_traits.h"

#ifdef MATRICE_SIMD_ARCH
#include <xmmintrin.h>
#include <immintrin.h>
#include <emmintrin.h>

#define MATRICE_SIMD_ARCH MATRICE_SIMD_AVX

MATRICE_ARCH_BEGIN
/// <summary>
/// \brief CONSTRAINT, define allowed types in SIMD.
/// </summary>
template<typename T>
concept packetable_scalar = is_any_of_v<T,
	float_t, double_t,
	int8_t, int16_t, int32_t, int64_t,
	uint8_t, uint16_t, uint32_t, uint64_t>;

/// <summary>
/// \brief Determine the allowed element size for SIMD vector.  
/// </summary>
template<packetable_scalar T>
struct packet_size {
	static constexpr uint8_t value =
#if MATRICE_SIMD_ARCH==MATRICE_SIMD_SSE      // __m128
	((1 << 7) >> 3) / sizeof(T);
#elif MATRICE_SIMD_ARCH==MATRICE_SIMD_AVX    // __m256
	((1 << 8) >> 3) / sizeof(T);
#elif MATRICE_SIMD_ARCH==MATRICE_SIMD_AVX512 // __m512
	((1 << 9) >> 3) / sizeof(T);
#else
		1;
#endif
};
/// <summary>
/// \brief Deduce the number of elements to be packed. 
/// </summary>
/// <typeparam name="T">Any scalar type supported in SIMD</typeparam>
template<packetable_scalar T>
inline constexpr auto packet_size_v = packet_size<T>::value;

/// <summary>
/// \brief Deduce vector data type in SIMD.
/// </summary>
template<packetable_scalar T, uint8_t N> struct packed_vector {
	static_assert(false, 
		"Parameter N must be 2, 4, 8, 16, 32, or 64 in packed_vector<T, N>.");
};
template<packetable_scalar T> struct packed_vector<T, 2> {
	using type = conditional_t<is_same_v<T, double>, __m128d, 
		conditional_t<is_any_of_v<T, int64_t, uint64_t>, __m128i, void>>;
};
template<packetable_scalar T> struct packed_vector<T, 4> {
	using type = conditional_t<is_same_v<T, float>, __m128, 
		conditional_t<is_same_v<T, double>, __m256d, 
		conditional_t<is_any_of_v<T, int64_t, uint64_t>, __m256i, 
		conditional_t<is_any_of_v<T, int32_t, uint32_t>, __m128i, void>>>>;
};
template<packetable_scalar T> struct packed_vector<T, 8> {
	using type = conditional_t<is_same_v<T, float>, __m256, 
		conditional_t<is_same_v<T, double>, __m512d, 
		conditional_t<is_any_of_v<T, int64_t, uint64_t>, __m512i,
		conditional_t<is_any_of_v<T, int32_t, uint32_t>, __m256i,  
		conditional_t<is_any_of_v<T, int16_t, uint16_t>, __m128i, void>>>>>;
};
template<packetable_scalar T> struct packed_vector<T, 16> {
	using type = conditional_t<is_same_v<T, float>, __m512, 
		conditional_t<is_any_of_v<T, int32_t, uint32_t>, __m512i,
		conditional_t<is_any_of_v<T, int16_t, uint16_t>, __m256i,  
		conditional_t<is_any_of_v<T, int8_t, uint8_t>, __m128i, void>>>>;
};
template<packetable_scalar T> struct packed_vector<T, 32> {
	using type = conditional_t<is_any_of_v<T, int16_t, uint16_t>, __m512i,
		conditional_t<is_any_of_v<T, int8_t, uint8_t>, __m256i, void>>;
};
template<packetable_scalar T> struct packed_vector<T, 64> {
	using type = conditional_t<is_any_of_v<T, int8_t, uint8_t>, __m512i, void>;
};
template<packetable_scalar T>
using packed_vector_t = packed_vector<T, packet_size_v<T>>::type;

/// <summary>
/// \brief TRAITS, for auto-deduce data type and size in SIMD.
/// </summary>
template<packetable_scalar T>
struct packet_traits : traits<T> {
	using value_type = remove_all_t<T>;
	using type = packed_vector_t<value_type>;
	static constexpr auto size = packet_size_v<value_type>;
};


MATRICE_ARCH_END
#endif // MATRICE_SIMD_ARCH
