/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include "../../private/_type_traits.h"
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

MATRICE_ARCH_BEGIN namespace detail{ namespace impl {
#ifndef MATRICE_HOST_ICEA
#define HOST_INL_CXPR_T MATRICE_HOST_FINL auto
#endif
#ifndef HOST_STATIC_INL_CXPR_T
#define HOST_STATIC_INL_CXPR_T static MATRICE_HOST_FINL constexpr auto 
#endif
#pragma region <!-- packet level operators (PLO) : data management & horizontal operation -->
template<typename T, int _Elems> struct simd_hop_base
{ 
	using value_t = T;
	using pointer = std::add_pointer_t<value_t>;
	using raw_type = conditional_t<value_t, _Elems>;

	template<size_t _N = _Elems - 1>  struct _Reduce_n {
		HOST_STATIC_INL_CXPR_T value(const pointer _Data[[_N]]) {
			return (_Reduce_n<_N - 1>::value(_Data) + _Data[_N]);
		}
	};
	template<> struct _Reduce_n<0> {
		HOST_STATIC_INL_CXPR_T value(const pointer _Data[[]]) {
			return (_Data[0]);
		}
	};
};
template<typename T, int _Elems> struct simd_hop : simd_hop_base<T, _Elems> {
	static_assert(true, "Oops! In simd_hop<T, _Elems>, T and/or _Elems may not be supported.");
};
template<> struct simd_hop<size_t, 4> : simd_hop_base<size_t, 4>
{
	HOST_INL_CXPR_T operator()(const value_t _Value) const noexcept
	{ return _mm256_set1_epi64x(_Value); }
	HOST_INL_CXPR_T operator()(const pointer _First) const noexcept
	{ return _mm256_load_si256((const raw_type*)_First); }
	HOST_INL_CXPR_T operator()(raw_type& _Packet) const noexcept
	{ return  _Packet.m256i_u64; }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet) const noexcept
	{ return  _Packet.m256i_u64; }
	HOST_INL_CXPR_T operator()(raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm256_store_si256((raw_type*)_Dst, _Packet); }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm256_store_si256((raw_type*)_Dst, _Packet); }
	HOST_INL_CXPR_T operator+ (const pointer _First) const noexcept
	{ return (_Reduce_n<>::value(_First)); }
};
template<> struct simd_hop<ptrdiff_t, 4> : simd_hop_base<ptrdiff_t, 4>
{
	HOST_INL_CXPR_T operator()(const value_t _Value) const noexcept
	{ return _mm256_set1_epi64x(_Value); }
	HOST_INL_CXPR_T operator()(const pointer _First) const noexcept
	{ return _mm256_load_si256((const raw_type*)_First); }
	HOST_INL_CXPR_T operator()(raw_type& _Packet) const noexcept
	{ return  _Packet.m256i_i64; }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet) const noexcept
	{ return  _Packet.m256i_i64; }
	HOST_INL_CXPR_T operator()(raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm256_store_si256((raw_type*)_Dst, _Packet); }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm256_store_si256((raw_type*)_Dst, _Packet); }
	HOST_INL_CXPR_T operator+ (const pointer _First) const noexcept
	{ return (_Reduce_n<>::value(_First)); }
};
template<> struct simd_hop<float, 4> : simd_hop_base<float, 4>
{
	HOST_INL_CXPR_T operator()(const value_t _Value) const noexcept
	{ return _mm_set_ps1(_Value); }
	HOST_INL_CXPR_T operator()(const pointer _First) const noexcept
	{ return _mm_load_ps(_First); }
	HOST_INL_CXPR_T operator()(raw_type& _Packet) const noexcept
	{ return  _Packet.m128_f32; }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet) const noexcept
	{ return  _Packet.m128_f32; }
	HOST_INL_CXPR_T operator()(raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm_store_ps(_Dst, _Packet); }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm_store_ps(_Dst, _Packet); }
	HOST_INL_CXPR_T operator+ (const pointer _First) const noexcept
	{ return (_Reduce_n<>::value(_First)); }
};
template<> struct simd_hop<float, 8> : simd_hop_base<float, 8>
{
	HOST_INL_CXPR_T operator()(const value_t _Value) const noexcept
	{ return _mm256_set1_ps(_Value); }
	HOST_INL_CXPR_T operator()(const pointer _First) const noexcept
	{ return _mm256_load_ps(_First); }
	HOST_INL_CXPR_T operator()(raw_type& _Packet) const noexcept
	{ return  _Packet.m256_f32; }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet) const noexcept
	{ return  _Packet.m256_f32; }
	HOST_INL_CXPR_T operator()(raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm256_store_ps(_Dst, _Packet); }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm256_store_ps(_Dst, _Packet); }
	HOST_INL_CXPR_T operator+ (const pointer _First) const noexcept
	{	return (_Reduce_n<>::value(_First)); }
};
template<> struct simd_hop<float, 16> : simd_hop_base<float, 16>
{
	HOST_INL_CXPR_T operator()(const value_t _Value) const noexcept
	{ return _mm512_set1_ps(_Value); }
	HOST_INL_CXPR_T operator()(const pointer _First) const noexcept
	{ return _mm512_load_ps(_First); }
	HOST_INL_CXPR_T operator()(raw_type& _Packet) const noexcept
	{ return  _Packet.m512_f32; }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet) const noexcept
	{ return  _Packet.m512_f32; }
	HOST_INL_CXPR_T operator()(raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm512_store_ps(_Dst, _Packet); }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm512_store_ps(_Dst, _Packet); }
	HOST_INL_CXPR_T operator+ (const pointer _First) const noexcept
	{	return (_Reduce_n<>::value(_First)); }
};
template<> struct simd_hop<double, 2> : simd_hop_base<double, 2>
{ 
	HOST_INL_CXPR_T operator()(const value_t _Value) const noexcept
	{ return _mm_set1_pd(_Value); }
	HOST_INL_CXPR_T operator()(const pointer _First) const noexcept
	{ return _mm_load_pd(_First); }
	HOST_INL_CXPR_T operator()(raw_type& _Packet) const noexcept
	{ return  _Packet.m128d_f64; }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet) const noexcept
	{ return  _Packet.m128d_f64; }
	HOST_INL_CXPR_T operator()(raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm_store_pd(_Dst, _Packet); }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm_store_pd(_Dst, _Packet); }
	HOST_INL_CXPR_T operator+ (const pointer _First) const noexcept
	{ return (_Reduce_n<>::value(_First)); }
};
template<> struct simd_hop<double, 4> : simd_hop_base<double, 4>
{
	HOST_INL_CXPR_T operator()(const value_t _Value) const noexcept
	{ return _mm256_set1_pd(_Value); }
	HOST_INL_CXPR_T operator()(const pointer _First) const noexcept
	{ return _mm256_load_pd(_First); }
	HOST_INL_CXPR_T operator()(raw_type& _Packet) const noexcept
	{ return  _Packet.m256d_f64; }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet) const noexcept
	{ return  _Packet.m256d_f64; }
	HOST_INL_CXPR_T operator()(raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm256_store_pd(_Dst, _Packet); }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm256_store_pd(_Dst, _Packet); }
	HOST_INL_CXPR_T operator+ (const pointer _First) const noexcept
	{ return (_Reduce_n<>::value(_First)); }
};
template<> struct simd_hop<double, 8> : simd_hop_base<double, 8>
{
	HOST_INL_CXPR_T operator()(const value_t _Value) const noexcept
	{ return _mm512_set1_pd(_Value); }
	HOST_INL_CXPR_T operator()(const pointer _First) const noexcept
	{ return _mm512_load_pd(_First); }
	HOST_INL_CXPR_T operator()(raw_type& _Packet) const noexcept
	{ return  _Packet.m512d_f64; }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet) const noexcept
	{ return  _Packet.m512d_f64; }
	HOST_INL_CXPR_T operator()(raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm512_store_pd(_Dst, _Packet); }
	HOST_INL_CXPR_T operator()(const raw_type& _Packet, pointer _Dst) const noexcept
	{ _mm512_store_pd(_Dst, _Packet); }
	HOST_INL_CXPR_T operator+ (const pointer _First) const noexcept
	{	return (_Reduce_n<>::value(_First)); }
};
#pragma endregion

#pragma region <!-- raw-type level operators (RLO) : vertical arithmetic operation -->
template<typename T, int _Elems> struct simd_vop { 
	enum { N = _Elems }; 
	using value_t = T; 
};

template<typename T, int _Elems, typename derived = simd_vop<T, _Elems>> struct simd_vop_base {
	enum { N = _Elems }; using value_t = T;
	using type = dgelom::simd::conditional_t<value_t, N>;
	template<typename _Op> HOST_STATIC_INL_CXPR_T _Binary(_Op _Op) noexcept { 
		return _Op(); }
	template<typename _Op> HOST_STATIC_INL_CXPR_T _Unary(_Op _Op) noexcept { 
		return _Op(); }
};

template<> struct simd_vop<size_t, 4> : public simd_vop_base<size_t, 4>
{
	using base_t = simd_vop_base<size_t, 4>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;

	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm256_add_epi64(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm256_sub_epi64(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm256_mul_epu32(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right) noexcept {
		//return base_t::_Binary([&]()->auto{return _mm256_div_epu32(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) noexcept {
		return base_t::_Unary([&]()->type {return _mm256_and_si256(_Right, _mm256_castsi128_si256(_mm_set1_epi32(~(1 << 31)))); });
	}
};
template<> struct simd_vop<float, 4> : public simd_vop_base<float, 4>
{
	using base_t = simd_vop_base<float, 4>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;

	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm_add_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm_sub_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm_mul_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm_div_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) noexcept {
		return  base_t::_Unary([&]()->type{return _mm_and_ps(_Right, _mm_castsi128_ps(_mm_set1_epi32(~(1 << 31)))); });
	}
};
template<> struct simd_vop<float, 8> : public simd_vop_base<float, 8>
{
	using base_t = simd_vop_base<float, 8>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;

	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm256_add_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm256_sub_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm256_mul_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm256_div_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) noexcept {
		return base_t::_Unary([&]()->type{ return _mm256_and_ps(_Right, _mm256_castsi256_ps(_mm256_set1_epi32(~(1 << 31)))); });
	}
};
template<> struct simd_vop<float, 16> : public simd_vop_base<float, 16>
{
	using base_t = simd_vop_base<float, 16>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;
	
	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm512_add_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm512_sub_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm512_mul_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm512_div_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) noexcept {
		return base_t::_Unary([&]()->type{return _mm512_castsi512_ps(_mm512_srli_epi64(_mm512_slli_epi64(_mm512_castps_si512(_Right), 1), 1)); });
	}
};
template<> struct simd_vop<double, 2> : public simd_vop_base<double, 2>
{
	using base_t = simd_vop_base<double, 2>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;
	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm_add_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm_sub_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm_mul_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm_div_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) noexcept {
		return base_t::_Unary([&]()->type{ return _mm_and_pd(_Right, _mm_castsi128_pd(_mm_setr_epi32(-1, 0x7FFFFFFF, -1, 0x7FFFFFFF)));});
	}
};
template<> struct simd_vop<double, 4> : public simd_vop_base<double, 4>
{
	using base_t = simd_vop_base<double, 4>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;

	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm256_add_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm256_sub_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm256_mul_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm256_div_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) noexcept {
		return base_t::_Unary([&]()->type{ return _mm256_and_pd(_Right, _mm256_castsi256_pd(_mm256_set1_epi32(~(1 << 31)))); });
	}
};
template<> struct simd_vop<double, 8> : public simd_vop_base<double, 8>
{
	using base_t = simd_vop_base<double, 8>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;

	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm512_add_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm512_sub_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm512_mul_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right) noexcept {
		return base_t::_Binary([&]()->auto{return _mm512_div_pd(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T const abs(const type& _Right) noexcept { return _Right; }
};
#pragma endregion

}} MATRICE_ARCH_END
#endif