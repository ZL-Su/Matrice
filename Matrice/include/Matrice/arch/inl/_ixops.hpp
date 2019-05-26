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
#include "../../private/_type_traits.h"
#ifdef __AVX__
#include "../_ixtraits.hpp"
#include "_ixop_impls.hpp"
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <wmmintrin.h>
#include <xmmintrin.h>

MATRICE_ARCH_BEGIN  namespace detail {

///<brief> SIMD operator definitions for application level. </brief>
template<
	typename T, 
	typename = std::enable_if_t<std::is_arithmetic_v<T>>>
struct Op_ MATRICE_NONHERITABLE {
	using value_t = T;
	template<size_t _Elems = 4> struct _base_type {
		enum {num_of_elem = _Elems};
		using type = conditional_t<value_t, num_of_elem>;
	};

	/**
	 * ... - vertical binary arithmetic operations.
	 */
	template<size_t _Elems> struct plus : _base_type<_Elems> {
		using type = typename _base_type<_Elems>::type;
		enum { size = _base_type<_Elems>::num_of_elem };
		HOST_INL_CXPR_T operator() (const type& _Left, const type& _Right) {
			return impl::simd_vop<value_t, size>::sum(_Left, _Right);
		}
	};
	template<size_t _Elems> struct minus {
		using type = typename _base_type<_Elems>::type;
		enum { size = _base_type<_Elems>::num_of_elem };
		HOST_INL_CXPR_T operator() (const type& _Left, const type& _Right) {
			return impl::simd_vop<value_t, size>::sub(_Left, _Right);
		}
	};
	template<size_t _Elems> struct multiplies {
		using type = typename _base_type<_Elems>::type;
		enum { size = _base_type<_Elems>::num_of_elem };
		HOST_INL_CXPR_T operator() (const type& _Left, const type& _Right) {
			return impl::simd_vop<value_t, size>::mul(_Left, _Right);
		}
	};
	template<size_t _Elems> struct divides {
		using type = typename _base_type<_Elems>::type;
		enum { size = _base_type<_Elems>::num_of_elem };
		HOST_INL_CXPR_T operator() (const type& _Left, const type& _Right) {
			return impl::simd_vop<value_t, size>::div(_Left, _Right);
		}
	};

	/**
	 * ... - vertical unary arithmetic operations.
	 */
	template<size_t _Elems> struct abs {
		using type = typename _base_type<_Elems>::type;
		enum { size = _base_type<_Elems>::num_of_elem };
		HOST_INL_CXPR_T operator() (const type& _Right) {
			return impl::simd_vop<value_t, size>::abs(_Right);
		}
	};

	/**
	 * adaptor - for data loading, storing and horizontal collection.
	 */
	template<size_t _Elems> using adaptor = impl::simd_hop<value_t, _Elems>;
};
template<typename _Fn, typename... _Args>
HOST_INL_CXPR_T _Transform_impl(const _Args&... _args) {
	//_Fn _Func;
	return _Fn()(_args...);
}
}
template<typename _Fn, typename... _Args>
HOST_INL_CXPR_T transform(const _Args&... _args) { 
	return (detail::_Transform_impl<_Fn>(_args.data()...)); 
}
MATRICE_ARCH_END

#endif