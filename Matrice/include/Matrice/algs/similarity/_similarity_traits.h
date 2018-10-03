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
template<typename _Ty, size_t _Opt> class _GaussNewton_impl_ic;

namespace internal {
// \size traits
template<size_t _Values> struct static_size { 
	enum { value = 
		(_Values & bilinear == bilinear) ? 3 :
		(_Values & bcspline == bcspline) ? 4 :
		(_Values & bqspline == bqspline) ? 6 :
		(_Values & bsspline == bsspline) ? 8 : _Values
	};
};

// \Gaussian-Newton solver traits
template<typename _Solver_type>
struct gn_solver_traits { 
	using type = typename _Solver_type::value_t;
	enum { options = _Solver_type::options };
};
template<typename _Ty, size_t _Options>
struct gn_solver_traits<_GaussNewton_impl_ic<_Ty, _Options>> {
	using type = _Ty;
	enum { options = _Options };
};
}}
MATRICE_ALGS_END
