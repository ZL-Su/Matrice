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
#include "../private/_expr_type_traits.h"
#include "../util/_macros.h"

#ifdef __AVX__
#include <pmmintrin.h>

#ifndef HOST_STATIC_INL_CXPR_T
#define HOST_STATIC_INL_CXPR_T static MATRICE_HOST_FINL constexpr auto 
#endif // !HOST_STATIC_INL_CXPR_T

MATRICE_ARCH_BEGIN
template<typename T, int _Elems> struct simd_op
{ enum { N = _Elems }; using value_t = T; };

template<typename T, int _Elems, typename derived = simd_op<T, _Elems>>
struct simd_op_base {
	enum { N = _Elems }; using value_t = T;
	using type = dgelom::simd::conditional_t<value_t, N>;
	template<typename _Op> static MATRICE_HOST_FINL auto _Binary(_Op _Op) { return _Op(); }
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
};
MATRICE_ARCH_END
#endif
