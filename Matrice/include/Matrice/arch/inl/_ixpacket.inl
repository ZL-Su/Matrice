
#pragma once
#include <pmmintrin.h>
#include "../ixpacket.h"
#ifndef HOST_STATIC_INL_CXPR_T
#define HOST_STATIC_INL_CXPR_T static MATRICE_HOST_FINL constexpr auto 
#endif // !HOST_STATIC_INL_CXPR_T

template<typename T, int _Elems> struct simd_op
{
	using value_t = T;
	enum { N = _Elems };
};

template<typename T, int _Elems, typename derived = simd_op<T, _Elems>>
struct simd_op_base
{
	enum { N = _Elems };
	using value_t = T;
	using type = dgelom::simd::conditional_t<value_t, N>;
	template<typename _Op> static 
	MATRICE_HOST_FINL auto _Binary(_Op _Op) { return _Op();}
};

template<> struct simd_op<float, 4> : public simd_op_base<float, 4>
{
	using base_t = simd_op_base<float, 4>;
	using base_t::value_t;
	using base_t::type;
	using base_t::N;
	HOST_STATIC_INL_CXPR_T sum(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm_hadd_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T sub(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm_hsub_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T mul(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm_mul_ps(_Left, _Right); });
	}
	HOST_STATIC_INL_CXPR_T div(const type& _Left, const type& _Right)
	{
		return base_t::_Binary([&]()->auto{return _mm_div_ps(_Left, _Right); });
	}
};