/*********************************************************************
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
***********************************************************************/
#include "../interpolation.h"
#pragma once

MATRICE_ALGS_BEGIN namespace detail {
// \forward declarations
template<
	typename _Data_type, 
	typename _Interpolation_tag, 
	std::size_t _Warp_fn_order     //can be 1 or 2
>  class _Invcomp_conv_impl;

namespace internal {
// \size traits
template<size_t _Values> struct static_size { 
	enum { value = 
		_Values == bilinear ? 3 : _Values == bcspline ? 4 :
		_Values == bqspline ? 6 : _Values == bsspline ? 8 : _Values
	};
};

// \iterative convolutional solver traits
template<typename _Solver_type> struct conv_solver_traits {};

template<typename _Ty, typename _Itag, std::size_t _O>
struct conv_solver_traits<_Invcomp_conv_impl<_Ty, _Itag, _O>> {
	using interp_category = _Itag;
	using value_type = _Ty;
	static const auto interp = INTERP|BSPLINE;
	static const auto order = _O;
};
}}
MATRICE_ALGS_END
