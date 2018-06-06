#pragma once
#include <initializer_list>
#include "../private/_expr_type_traits.h"
#include "../util/_macros.h"
#include "../util/genalgs.h"
#ifdef __AVX__
#include <mmintrin.h>    //_m64
#include <emmintrin.h>   //_m128
#include <immintrin.h>   //_m256

#define MATRICE_ARCH_BEGIN MATRICE_NAMESPACE_BEGIN_ namespace simd {
#define MATRICE_ARCH_END   } _MATRICE_NAMESPACE_END

MATRICE_ARCH_BEGIN
#pragma region <!-- simd traits -->
template<typename T, int _Elems> struct conditional {};
template<typename T> struct conditional<T, 4>
{
	using type = typename dgelom::conditional<dgelom::is_float32<T>::value, __m128, typename dgelom::conditional<dgelom::is_float64<T>::value, __m256d, void>::type>::type;
};
template<typename T> struct conditional<T, 8>
{
	using type = typename dgelom::conditional<dgelom::is_float32<T>::value, __m256, typename dgelom::conditional<dgelom::is_float64<T>::value, __m512d, void>::type>::type;
};
template<typename T, int _Elems> using conditional_t = typename conditional<T, _Elems>::type;

template<typename T, int _Elems> struct simd_traits
{
	struct is_sfp8 { enum { value = is_float32<T>::value && _Elems == 8 }; };
	struct is_dfp8 { enum { value = is_float64<T>::value && _Elems == 8 }; };
};
#pragma endregion

#pragma region <!-- packet setters -->
template<typename T, int _Elems> struct set_packet { using value_t = T; };
template<> struct set_packet<float,  1<<1<<1>
{
	using value_t = float; using pointer = value_t*;
	MATRICE_HOST_INL constexpr auto operator()(const value_t _Value) const
	{ return _mm_set_ps1(_Value); }
	MATRICE_HOST_INL constexpr auto operator()(pointer _First) const
	{ return _mm_set_ps(*(_First), *(_First++), *(_First++), *(_First++)); }
};
template<> struct set_packet<float,  1<<1<<1<<1>
{
	using value_t = float; using pointer = value_t*;
	MATRICE_HOST_INL constexpr auto operator()(const value_t _Value) const
	{ return _mm256_set1_ps(_Value); }
	MATRICE_HOST_INL constexpr auto operator()(pointer _First) const
	{ return _mm256_set_ps(*(_First), *(_First++), *(_First++), *(_First++), *(_First++), *(_First++), *(_First++), *(_First++)); }
};
template<> struct set_packet<float,  1<<1<<1<<1<<1>
{
	using value_t = float; using pointer = value_t*;
	MATRICE_HOST_INL constexpr auto operator()(const value_t _Value) const
	{ return _mm256_set1_ps(_Value); }
	MATRICE_HOST_INL constexpr auto operator()(pointer _First) const
	{ return _mm256_set_ps(*(_First), *(_First++), *(_First++), *(_First++), *(_First++), *(_First++), *(_First++), *(_First++)); }
};
template<> struct set_packet<double, 1<<1>
{
	using value_t = float; using pointer = value_t*;
	MATRICE_HOST_INL constexpr auto operator()(const value_t _Value) const
	{ return _mm_set1_pd(_Value); }
	MATRICE_HOST_INL constexpr auto operator()(pointer _First) const
	{ return _mm_set_pd(*(_First), *(_First++)); }
};
template<> struct set_packet<double, 1<<1<<1>
{
	using value_t = float; using pointer = value_t*;
	MATRICE_HOST_INL constexpr auto operator()(const value_t _Value) const
	{ return _mm256_set1_pd(_Value); }
	MATRICE_HOST_INL constexpr auto operator()(pointer _First) const
	{ return _mm256_set_pd(*(_First), *(_First++), *(_First++), *(_First++)); }
};
template<> struct set_packet<double, 1<<1<<1<<1>
{
	using value_t = float; using pointer = value_t*;
	MATRICE_HOST_INL constexpr auto operator()(const value_t _Value) const
	{ return _mm512_set1_pd(_Value); }
	MATRICE_HOST_INL constexpr auto operator()(pointer _First) const
	{ return _mm512_set_pd(*(_First), *(_First++), *(_First++), *(_First++), *(_First++), *(_First++), *(_First++), *(_First++)); }
};
#pragma endregion

template<typename T, int _Elems> class simd_base_
{
	set_packet<T, _Elems> _Set_op;
public:
	using value_t    = T;
	using pointer    = value_t*;
	using const_pointer = const pointer;
	using internal_t = conditional_t<value_t, _Elems>;
	using initlist_t = std::initializer_list<value_t>;

protected:
	internal_t m_data;
	template<typename... _Args> MATRICE_HOST_FINL 
	void _Set(_Args... _args) { m_data = _Set_op(_args...); }
};

//template<typename T, int _Elems, typename Packet = conditional_t<T, _Elems>>
MATRICE_HOST_INL auto operator+ (const conditional_t<float, 8>& _Left, const conditional_t<float, 8>& _Right)
{
	//if constexpr (simd_traits<T, _Elems>::is_sfp8::value)
		return _mm256_add_ps(_Left, _Right);
	//if constexpr (simd_traits<T, _Elems>::is_dfp8::value)
		//return _mm512_add_pd(_Left, _Right);
}
template<typename T, int _Elems, typename Packet = conditional_t<T, _Elems>>
MATRICE_HOST_INL auto operator- (const Packet& _Left, const Packet& _Right)
{
	if constexpr (simd_traits<T, _Elems>::is_sfp8::value)
		return _mm256_sub_ps(_Left, _Right);
	if constexpr (simd_traits<T, _Elems>::is_dfp8::value)
		return _mm512_sub_pd(_Left, _Right);
}
template<typename T, int _Elems, typename Packet = conditional_t<T, _Elems>>
MATRICE_HOST_INL auto operator* (const Packet& _Left, const Packet& _Right)
{
	if constexpr (simd_traits<T, _Elems>::is_sfp8::value)
		return _mm256_mul_ps(_Left, _Right);
	if constexpr (simd_traits<T, _Elems>::is_dfp8::value)
		return _mm512_mul_pd(_Left, _Right);
}
template<typename T, int _Elems, typename Packet = conditional_t<T, _Elems>>
MATRICE_HOST_INL auto operator/ (const Packet& _Left, const Packet& _Right)
{
	if constexpr (simd_traits<T, _Elems>::is_sfp8::value)
		return _mm256_div_ps(_Left, _Right);
	if constexpr (simd_traits<T, _Elems>::is_dfp8::value)
		return _mm512_div_pd(_Left, _Right);
}
template<typename T, int _Elems> MATRICE_HOST_INL
auto reduce(conditional_t<T, _Elems> packet, std::enable_if_t<is_float32<T>::value>* = 0)
{
	return dgelom::reduce(&packet[0], &packet[0]+_Elems);
}
MATRICE_ARCH_END
#endif
