/*  *************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
*	*************************************************************************/
#pragma once
#include <type_traits>
#include "../util/_type_defs.h"
#include "_memory.h"

namespace dgelom {
template<typename T> struct traits {};
template<typename T> struct traits<const T> : traits<T> { using type = T; };
struct MatrixExpr {};
template<typename Derived, typename Expression> struct DenseExprBase {};
template<typename Derived> struct DenseExprBase<Derived, MatrixExpr> {};
template<typename _Ty> struct remove_reference { using type = typename std::remove_reference<_Ty>::type; };
template<typename _Ty> struct type_bytes { enum { value = sizeof(_Ty) }; };
template<bool _Test, typename T1, typename T2> struct conditonal {};
template<typename T1, typename T2> struct conditonal<true, T1, T2> { using type = T1; };
template<typename T1, typename T2> struct conditonal<false, T1, T2> { using type = T2; };
template<bool _Test, typename T1, typename T2> struct conditional {};
template<typename T1, typename T2> struct conditional<true, T1, T2> { using type = T1; };
template<typename T1, typename T2> struct conditional<false, T1, T2> { using type = T2; };
template<bool _Test, typename T1, typename T2> using conditional_t = typename conditional<_Test, T1, T2>::type;
template<int _Opt> struct is_expression { enum { value = _Opt & expr == expr ? true : false }; };
template<int _Val> struct is_zero { enum { value = _Val == 0 ? true : false }; };
template<int _R, int _C> struct is_static {enum {value = _R > 0 && _C >0 ? true : false}; };
template<typename T> struct is_common_int64 { enum { value = std::is_integral_v<T> && sizeof(T) == 8 }; };
template<typename T> struct is_int64 { enum { value = std::is_signed_v<T> && std::is_integral_v<T> && sizeof(T) == 8 }; };
template<typename T> struct is_uint64 { enum { value = std::is_unsigned_v<T> && std::is_integral_v<T> && sizeof(T) == 8 }; };
template<typename T> struct is_float32 { enum { value = std::is_floating_point<T>::value && (sizeof(T) == 4) }; };
template<typename T> struct is_float64 { enum { value = std::is_floating_point<T>::value && (sizeof(T) == 8) }; };
template<int _M, int _N> struct allocator_option{ 
	enum {
		value = _M >0  && _N>0   ? LINEAR + COPY :  // stack allocator
		        _M==0  && _N==-1 ? LINEAR        :  // linear device allocator
		        _M==-1 && _N==-1 ? PITCHED       :  // pitched device allocator
#ifdef __CXX11_SHARED__
										 LINEAR + SHARED  // smart heap or global allocator
#else
		                         LINEAR + COPY    // deep heap or global allocator
#endif      
	}; 
};

template<int _Rows = 0, int _Cols = 0> struct compile_time_size {
	enum { val_1 = 0x0001, val_2 = 0x0002, val_3 = 0x0003, val_4 = 0x0004 };
	enum {
		RunTimeDeduceInHost = 0,
		RunTimeDeduceInDevi = -1,
		CompileTimeRows = _Rows,
		CompileTimeCols = _Cols,
	};
	static const int RunTimeDeducedInHost = 0, RunTimeDeducedInDevice = -1;
};
}