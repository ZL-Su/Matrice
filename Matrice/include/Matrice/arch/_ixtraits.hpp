#pragma once
#include "../private/_expr_type_traits.h"
#include "../util/_macros.h"
#ifdef __AVX__
#include <mmintrin.h>    //_m64
#include <emmintrin.h>   //_m128
#include <immintrin.h>   //_m256
#include <zmmintrin.h>   //_m512
MATRICE_ARCH_BEGIN
#pragma region <!-- simd traits -->
template<typename T, int _Elems> struct conditional {};
template<typename T> struct conditional<T, 2>
{
	using type = dgelom::conditional_t<is_common_int64<T>::value, __m128i, dgelom::conditional_t<is_float32<T>::value, __m64, dgelom::conditional_t<is_float64<T>::value, __m128d, void>>>;
};
template<typename T> struct conditional<T, 4>
{
	using type = dgelom::conditional_t<is_common_int64<T>::value, __m256i, dgelom::conditional_t<is_float32<T>::value, __m128, dgelom::conditional_t<is_float64<T>::value, __m256d, void>>>;
};
template<typename T> struct conditional<T, 8>
{
	using type = dgelom::conditional_t<is_common_int64<T>::value, __m512i, dgelom::conditional_t<is_float32<T>::value, __m256, dgelom::conditional_t<is_float64<T>::value, __m512d, void>>>;
};
template<typename T> struct conditional<T, 16>
{
	using type = typename dgelom::conditional<dgelom::is_float32<T>::value, __m512, typename dgelom::conditional<dgelom::is_float64<T>::value, __m512d, void>::type>::type;
};
template<typename T, int _Elems> using conditional_t = typename conditional<T, _Elems>::type;

template<typename T, int _Elems> struct simd_traits
{
	struct is_sfp8 { enum { value = is_float32<T>::value && _Elems == 8 }; };
	struct is_dfp8 { enum { value = is_float64<T>::value && _Elems == 8 }; };
};
#pragma endregion
MATRICE_ARCH_END
#endif // __AVX__