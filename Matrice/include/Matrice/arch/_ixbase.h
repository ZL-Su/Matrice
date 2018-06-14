/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
#include <initializer_list>
#include "../private/_expr_type_traits.h"
#include "../util/_macros.h"
#include "../util/genalgs.h"
#ifdef __AVX__
#include <mmintrin.h>    //_m64
#include <emmintrin.h>   //_m128
#include <immintrin.h>   //_m256
#include <zmmintrin.h>   //_m512

#define MATRICE_ARCH_BEGIN MATRICE_NAMESPACE_BEGIN_ namespace simd {
#define MATRICE_ARCH_END   } _MATRICE_NAMESPACE_END

#define MATRICE_HOST_ICEA MATRICE_HOST_INL constexpr auto

MATRICE_ARCH_BEGIN
#pragma region <!-- simd traits -->
template<typename T, int _Elems> struct conditional {};
template<typename T> struct conditional<T, 2>
{
	using type = typename dgelom::conditional<dgelom::is_float32<T>::value, __m64, typename dgelom::conditional<dgelom::is_float64<T>::value, __m128d, void>::type>::type;
};
template<typename T> struct conditional<T, 4>
{
	using type = typename dgelom::conditional<dgelom::is_float32<T>::value, __m128, typename dgelom::conditional<dgelom::is_float64<T>::value, __m256d, void>::type>::type;
};
template<typename T> struct conditional<T, 8>
{
	using type = typename dgelom::conditional<dgelom::is_float32<T>::value, __m256, typename dgelom::conditional<dgelom::is_float64<T>::value, __m512d, void>::type>::type;
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

#pragma region <!-- packet setters -->
template<typename T, int _Elems> struct packet_op { using value_t = T; }; 
template<> struct packet_op<float, 4>
{
	using value_t = float; using pointer = value_t*;
	using raw_packt_t = conditional_t<float, 4>;
	MATRICE_HOST_ICEA operator()(const value_t _Value) const
	{ return _mm_set_ps1(_Value); }
	MATRICE_HOST_ICEA operator()(const pointer _First) const
	{ return _mm_load_ps(_First); }
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet) const
	{ return  _Packet.m128_f32; }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet) const { return  _Packet.m128_f32; }
};
template<> struct packet_op<float, 8>
{
	using value_t = float; using pointer = value_t*;
	using raw_packt_t = conditional_t<float, 8>;
	MATRICE_HOST_ICEA operator()(const value_t _Value) const
	{ return _mm256_set1_ps(_Value); }
	MATRICE_HOST_ICEA operator()(const pointer _First) const
	{ return _mm256_load_ps(_First); }
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet) const
	{ return  _Packet.m256_f32; }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet) const
	{ return  _Packet.m256_f32; }
};
template<> struct packet_op<float, 16>
{
	using value_t = float; using pointer = value_t*;
	using raw_packt_t = conditional_t<float, 16>;
	MATRICE_HOST_ICEA operator()(const value_t _Value) const
	{ return _mm512_set1_ps(_Value); }
	MATRICE_HOST_ICEA operator()(const pointer _First) const
	{ return _mm512_load_ps(_First); }
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet) const
	{ return  _Packet.m512_f32; }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet) const
	{ return  _Packet.m512_f32; }
};
template<> struct packet_op<double, 2>
{
	using value_t = double; using pointer = value_t*;
	using raw_packt_t = conditional_t<double, 2>;
	MATRICE_HOST_ICEA operator()(const value_t _Value) const
	{ return _mm_set1_pd(_Value); }
	MATRICE_HOST_ICEA operator()(const pointer _First) const
	{ return _mm_load_pd(_First); }
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet) const
	{ return  _Packet.m128d_f64; }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet) const
	{ return  _Packet.m128d_f64; }
};
template<> struct packet_op<double, 4>
{
	using value_t = double; using pointer = value_t*;
	using raw_packt_t = conditional_t<double, 4>;
	MATRICE_HOST_ICEA operator()(const value_t _Value) const
	{ return _mm256_set1_pd(_Value); }
	MATRICE_HOST_ICEA operator()(const pointer _First) const
	{ return _mm256_load_pd(_First); }
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet) const
	{ return  _Packet.m256d_f64; }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet) const
	{ return  _Packet.m256d_f64; }
};
template<> struct packet_op<double, 8>
{
	using value_t = double; using pointer = value_t*;
	using raw_packt_t = conditional_t<double, 8>;
	MATRICE_HOST_ICEA operator()(const value_t _Value) const
	{ return _mm512_set1_pd(_Value); }
	MATRICE_HOST_ICEA operator()(const pointer _First) const
	{ return _mm512_load_pd(_First); }
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet) const
	{ return  _Packet.m512d_f64; }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet) const
	{ return  _Packet.m512d_f64; }
};
#pragma endregion

template<typename T, int _Elems> class simd_base_
{
	packet_op<T, _Elems> _op;
public:
	using value_t    = T;
	using const_value_t = const value_t;
	using pointer    = value_t*;
	using const_pointer = const pointer;
	using internal_t = conditional_t<value_t, _Elems>;
	using initlist_t = std::initializer_list<value_t>;
	MATRICE_HOST_FINL simd_base_(const_value_t _arg) noexcept { m_data = _op(_arg); }
	MATRICE_HOST_FINL simd_base_(const_pointer _arg) noexcept { m_data = _op(_arg); }

	MATRICE_HOST_FINL constexpr auto& operator[](size_t i) { return _op(m_data)[i]; }
	MATRICE_HOST_FINL constexpr const auto& operator[](size_t i) const { return _op(m_data)[i]; }
	MATRICE_HOST_FINL auto begin() { return _op(m_data); }
	MATRICE_HOST_FINL const auto begin() const { return _op(m_data); }
	MATRICE_HOST_FINL auto end() { return _op(m_data)+_Elems; }
	MATRICE_HOST_FINL const auto end() const { return _op(m_data) + _Elems; }
	MATRICE_HOST_INL const auto reduce() const { return dgelom::reduce<value_t>(begin(), end()); }

protected:
	template<typename... _Args> MATRICE_HOST_FINL constexpr
	auto _Op(_Args... _args)  { return _op(_args...); }
	template<typename... _Args> MATRICE_HOST_FINL constexpr
	auto _Op(_Args&... _args) { return _op(_args...); }
	internal_t m_data;
};

template<typename T, int _Elems> MATRICE_HOST_INL 
T reduce(const simd_base_<T, _Elems>& _Packed)
{
	return dgelom::reduce(_Packed.begin(), _Packed.end());
}

MATRICE_ARCH_END
#endif
