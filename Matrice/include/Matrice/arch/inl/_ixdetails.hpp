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
#include "../../private/_expr_type_traits.h"
#include "../../util/_macros.h"

#ifdef __AVX__
#include <mmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <wmmintrin.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <zmmintrin.h>
#include <xmmintrin.h>

MATRICE_ARCH_BEGIN namespace details{ namespace impl {
#ifndef MATRICE_HOST_ICEA
#define MATRICE_HOST_ICEA MATRICE_HOST_FINL constexpr auto
#endif
#ifndef HOST_STATIC_INL_CXPR_T
#define HOST_STATIC_INL_CXPR_T static MATRICE_HOST_FINL constexpr auto 
#endif
#pragma region <!-- packet level operators (PLO) : data management & horizontal operation -->
template<typename T, int _Elems> struct packet_op_base 
{ 
	using value_t = T;
	using pointer = value_t*;
	using raw_type = conditional_t<value_t, _Elems>;
};
template<typename T, int _Elems> 
struct packet_op : packet_op_base<T, _Elems>
{
	static_assert(true, "Oops, unsupported type.");
};
template<> struct packet_op<size_t, 4> : packet_op_base<size_t, 4>
{
	MATRICE_HOST_ICEA operator()(const value_t _Value) const
	{ return _mm256_set1_epi64x(_Value); }
	MATRICE_HOST_ICEA operator()(const pointer _First) const
	{ return _mm256_load_si256((const raw_type*)_First); }
	MATRICE_HOST_ICEA operator()(raw_type& _Packet) const
	{ return  _Packet.m256i_u64; }
	MATRICE_HOST_ICEA operator()(const raw_type& _Packet) const
	{ return  _Packet.m256i_u64; }
	MATRICE_HOST_ICEA operator()(raw_type& _Packet, pointer _Dst) const
	{ _mm256_store_si256((raw_type*)_Dst, _Packet); }
	MATRICE_HOST_ICEA operator()(const raw_type& _Packet, pointer _Dst) const
	{ _mm256_store_si256((raw_type*)_Dst, _Packet); }
	MATRICE_HOST_ICEA operator+ (const pointer _First) const
	{ return (_First[0] + _First[1] + _First[2] + _First[3]); }
};
template<> struct packet_op<ptrdiff_t, 4> : packet_op_base<ptrdiff_t, 4>
{
	MATRICE_HOST_ICEA operator()(const value_t _Value) const
	{ return _mm256_set1_epi64x(_Value); }
	MATRICE_HOST_ICEA operator()(const pointer _First) const
	{ return _mm256_load_si256((const raw_type*)_First); }
	MATRICE_HOST_ICEA operator()(raw_type& _Packet) const
	{ return  _Packet.m256i_i64; }
	MATRICE_HOST_ICEA operator()(const raw_type& _Packet) const
	{ return  _Packet.m256i_i64; }
	MATRICE_HOST_ICEA operator()(raw_type& _Packet, pointer _Dst) const
	{ _mm256_store_si256((raw_type*)_Dst, _Packet); }
	MATRICE_HOST_ICEA operator()(const raw_type& _Packet, pointer _Dst) const
	{ _mm256_store_si256((raw_type*)_Dst, _Packet); }
	MATRICE_HOST_ICEA operator+ (const pointer _First) const
	{ return (_First[0] + _First[1] + _First[2] + _First[3]); }
};
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
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet) const 
	{ return  _Packet.m128_f32; }
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet, pointer _Dst) const
	{ _mm_store_ps(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet, pointer _Dst) const
	{ _mm_store_ps(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator+ (const pointer _First) const
	{ return (_First[0] + _First[1] + _First[2] + _First[3]); }
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
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet, pointer _Dst) const
	{ _mm256_store_ps(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet, pointer _Dst) const
	{ _mm256_store_ps(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator+ (const pointer _First) const
	{	return (_First[0] + _First[1] + _First[2] + _First[3] + _First[4] + _First[5] + _First[6] + _First[7]);}
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
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet, pointer _Dst) const
	{ _mm512_store_ps(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet, pointer _Dst) const
	{ _mm512_store_ps(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator+ (const pointer _First) const
	{	return (_First[0] + _First[1] + _First[2] + _First[3] + _First[4] + _First[5] + _First[6] + _First[7]+ _First[8] + _First[9] + _First[10] + _First[11] + _First[12] + _First[13] + _First[14] + _First[15]);}
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
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet, pointer _Dst) const
	{ _mm_store_pd(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet, pointer _Dst) const
	{ _mm_store_pd(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator+ (const pointer _First) const
	{ return (_First[0] + _First[1]); }
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
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet, pointer _Dst) const
	{ _mm256_store_pd(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet, pointer _Dst) const
	{ _mm256_store_pd(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator+ (const pointer _First) const
	{ return (_First[0] + _First[1] + _First[2] + _First[3]); }
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
	MATRICE_HOST_ICEA operator()(raw_packt_t& _Packet, pointer _Dst) const
	{ _mm512_store_pd(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator()(const raw_packt_t& _Packet, pointer _Dst) const
	{ _mm512_store_pd(_Dst, _Packet); }
	MATRICE_HOST_ICEA operator+ (const pointer _First) const
	{	return (_First[0] + _First[1] + _First[2] + _First[3] + _First[4] + _First[5] + _First[6] + _First[7]);}
};
#pragma endregion

#pragma region <!-- raw-type level operators (RLO) : vertical arithmetic operation -->
template<typename T, int _Elems> struct simd_op
{ enum { N = _Elems }; using value_t = T; };

template<typename T, int _Elems, typename derived = simd_op<T, _Elems>>
struct simd_op_base {
	enum { N = _Elems }; using value_t = T;
	using type = dgelom::simd::conditional_t<value_t, N>;
	template<typename _Op> static MATRICE_HOST_FINL auto _Binary(_Op _Op) { return _Op(); }
	template<typename _Op> static MATRICE_HOST_FINL auto _Unary(_Op _Op) { return _Op(); }
};

template<> struct simd_op<size_t, 4> : public simd_op_base<size_t, 4>
{
	using base_t = simd_op_base<size_t, 4>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;
	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm256_add_epi64(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm256_sub_epi64(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm256_mul_epu32(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right)
	{
		//return base_t::_Binary([&]()->auto{return _mm256_div_epu64(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right)
	{
		return base_t::_Unary([&]()->type {return _mm256_and_si256(_Right, _mm256_castsi128_si256(_mm_set1_epi32(~(1 << 31)))); });
	}
};

template<> struct simd_op<float, 4> : public simd_op_base<float, 4>
{
	using base_t = simd_op_base<float, 4>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;
	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm_add_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm_sub_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm_mul_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right)
	{  
		return base_t::_Binary([&]()->auto{return _mm_div_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) 
	{
		return  base_t::_Unary([&]()->type{return _mm_and_ps(_Right, _mm_castsi128_ps(_mm_set1_epi32(~(1 << 31)))); });
	}
};
template<> struct simd_op<float, 8> : public simd_op_base<float, 8>
{
	using base_t = simd_op_base<float, 8>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;

	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm256_add_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm256_sub_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm256_mul_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm256_div_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right)
	{
		return base_t::_Unary([&]()->type{ return _mm256_and_ps(_Right, _mm256_castsi256_ps(_mm256_set1_epi32(~(1 << 31)))); });
	}
};
template<> struct simd_op<float, 16> : public simd_op_base<float, 16>
{
	using base_t = simd_op_base<float, 16>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;
	
	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm512_add_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm512_sub_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm512_mul_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm512_div_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) 
	{
		return base_t::_Unary([&]()->type{return _mm512_castsi512_ps(_mm512_srli_epi64(_mm512_slli_epi64(_mm512_castps_si512(_Right), 1), 1)); });
	}
};
template<> struct simd_op<double, 2> : public simd_op_base<double, 2>
{
	using base_t = simd_op_base<double, 2>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;
	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm_add_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm_sub_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm_mul_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm_div_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) 
	{ 
		return base_t::_Unary([&]()->type{ return _mm_and_pd(_Right, _mm_castsi128_pd(_mm_setr_epi32(-1, 0x7FFFFFFF, -1, 0x7FFFFFFF)));});
	}
};
template<> struct simd_op<double, 4> : public simd_op_base<double, 4>
{
	using base_t = simd_op_base<double, 4>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;

	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm256_add_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm256_sub_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm256_mul_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm256_div_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) 
	{
		return base_t::_Unary([&]()->type{ return _mm256_and_pd(_Right, _mm256_castsi256_pd(_mm256_set1_epi32(~(1 << 31)))); });
	}
};
template<> struct simd_op<double, 8> : public simd_op_base<double, 8>
{
	using base_t = simd_op_base<double, 8>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;

	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm512_add_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm512_sub_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm512_mul_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm512_div_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) { return _Right; }
};
#pragma endregion

}} MATRICE_ARCH_END
#endif