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

#include "../util/_macros.h"

DGE_MATRICE_BEGIN

// \default vector unit size for SIMD support
constexpr std::size_t packet_size_v =
#ifdef __AVX__
	4
#elif  __AVX2__
	8
#endif
;

template<int _Rows = 0, int _Cols = 0> struct compile_time_size {
	enum { val_1 = 0x0001, val_2 = 0x0002, val_3 = 0x0003, val_4 = 0x0004 };
	enum { CompileTimeRows = _Rows,  CompileTimeCols = _Cols,
			 RunTimeDeducedOnHost = 0, RunTimeDeducedOnDevice = -1 };
	static constexpr auto _1 = 0x0001;
	static constexpr auto _2 = 0x0002;
	static constexpr auto _3 = 0x0003;
	static constexpr auto _4 = 0x0004;
	static constexpr auto _5 = 0x0005;
	static constexpr auto _6 = 0x0006;
};
// \compile-time size of row value
template<int _M, int _N> MATRICE_GLOBAL_INL constexpr int ct_size_rv = compile_time_size<_M, _N>::CompileTimeRows;
// \compile-time size of column value
template<int _M, int _N> MATRICE_GLOBAL_INL constexpr int ct_size_cv = compile_time_size<_M, _N>::CompileTimeCols;
// \statically deduced runtime-size on host
template<int _M, int _N> MATRICE_GLOBAL_INL constexpr int rs_host_v = compile_time_size<_M, _N>::RunTimeDeducedOnHost;
// \statically deduced runtime-size on device
template<int _M, int _N> MATRICE_GLOBAL_INL constexpr int rs_device_v = compile_time_size<_M, _N>::RunTimeDeducedOnDevice;


template<bool _Test, int _N1, int _N2> struct conditional_size {};
template<int _N1, int _N2 > struct conditional_size<std::true_type::value, _N1, _N2> { enum { value = _N1 }; };
template<int _N1, int _N2 > struct conditional_size<std::false_type::value, _N1, _N2> { enum { value = _N2 }; };
template<bool _Test, int _N1, int _N2> MATRICE_GLOBAL_INL constexpr int conditional_size_v = conditional_size<_Test, _N1, _N2>::value;


template<int _N1, int _N2> struct max_integer { enum { value = conditional_size_v<(_N1 > _N2), _N1, _N2> }; };
template<int _N1, int _N2> MATRICE_GLOBAL_INL constexpr int max_integer_v = max_integer<_N1, _N2>::value;
template<int _N1, int _N2> struct min_integer { enum { value = conditional_size_v<(_N1 > _N2), _N2, _N1> }; };
template<int _N1, int _N2> MATRICE_GLOBAL_INL constexpr int min_integer_v = min_integer<_N1, _N2>::value;

DGE_MATRICE_END
