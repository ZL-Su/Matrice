
#pragma once
#include <pmmintrin.h>
#include "../ixpacket.h"
#ifndef HOST_STATIC_INL_CXPR_T
#define HOST_STATIC_INL_CXPR_T static MATRICE_HOST_FINL constexpr auto 
#endif // !HOST_STATIC_INL_CXPR_T

template<typename derived> struct simd_op_base
{
	using type = dgelom::simd::conditional_t<typename derived::value_t, derived::N>;
	template<typename _Op> static MATRICE_HOST_FINL auto _Binary(_Op _Op/*, const type& _Left, const type& _Right*/)
	{
		return _Op(/*_Left, _Right*/);
	}
};
template<typename T, int _Elems> struct simd_op 
{ 
	using value_t = T; 
	enum { N = _Elems };
};
template<> struct simd_op<float, 4> : public simd_op_base<simd_op<float, 4>>
{
	enum { N = 4 };
	using value_t = float;
	using base_t  = simd_op_base<simd_op<value_t, N>>;
	using type    = typename base_t::type;
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