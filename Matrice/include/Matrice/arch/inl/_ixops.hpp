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
#ifdef __AVX__
#include "../_ixtraits.hpp"
#include "_ixdetails.hpp"
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <wmmintrin.h>
#include <xmmintrin.h>

MATRICE_ARCH_BEGIN 
namespace details {
#define matrice_inl_cxauto MATRICE_HOST_FINL constexpr auto

template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
struct Op_ MATRICE_NONHERITABLE
{
	using value_t = T;
	template<size_t _Elems> struct _base_type
	{
		enum {num_of_elem = _Elems};
		using type = conditional_t<value_t, num_of_elem>;
	};
	template<size_t _Elems> struct plus
	{
		enum { num_of_elem = _Elems };
		using type = conditional_t<value_t, num_of_elem>;
		matrice_inl_cxauto operator() (const type& _Left, const type& _Right)
		{
			return impl::simd_op<value_t, num_of_elem>::sum(_Left, _Right);
		}
	};
	template<size_t _Elems> struct minus
	{
		enum { num_of_elem = _Elems };
		using type = conditional_t<value_t, num_of_elem>;
		matrice_inl_cxauto operator() (const type& _Left, const type& _Right)
		{
			return impl::simd_op<value_t, num_of_elem>::sub(_Left, _Right);
		}
	};
	template<size_t _Elems> struct multiplies
	{
		enum { num_of_elem = _Elems };
		using type = conditional_t<value_t, num_of_elem>;
		matrice_inl_cxauto operator() (const type& _Left, const type& _Right)
		{
			return impl::simd_op<value_t, num_of_elem>::mul(_Left, _Right);
		}
	};
	template<size_t _Elems> struct divides
	{
		enum { num_of_elem = _Elems };
		using type = conditional_t<value_t, num_of_elem>;
		matrice_inl_cxauto operator() (const type& _Left, const type& _Right)
		{
			return impl::simd_op<value_t, num_of_elem>::div(_Left, _Right);
		}
	};
	template<size_t _Elems> struct abs
	{
		enum { num_of_elem = _Elems };
		using type = conditional_t<value_t, num_of_elem>;
		matrice_inl_cxauto operator() (const type& _Right)
		{
			return impl::simd_op<value_t, num_of_elem>::abs(_Right);
		}
	};
};
template<typename _Fn, typename... _Args>
matrice_inl_cxauto _Transform_impl(_Fn _Func, const _Args&... _args) { return _Func(_args...); }
}
template<typename _Fn, typename... _Args>
matrice_inl_cxauto transform(_Fn _Func, const _Args&... _args) { return (details::_Transform_impl(_Func, _args.data()...)); }
MATRICE_ARCH_END

#endif