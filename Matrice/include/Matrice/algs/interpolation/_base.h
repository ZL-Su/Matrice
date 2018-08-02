/***********************************************************************
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
#pragma once

#include "../../core/matrix.h"
#include "../../core/vector.h"
#include "../../core/solver.h"

MATRICE_ALGS_BEGIN
enum {
	INTERP     = 18,

	BSPLINE    = 2307,
	BILINEAR   = 2152,
	BICUBIC    = 2156,
	BIQUINTIC  = 2160,
	BISEPTIC   = 2164
};

template<typename _Ty, size_t _Options> class BilinearInterp;
template<typename _Ty, size_t _Options> class BicubicSplineInterp;
template<typename _Ty, size_t _Options> class BiquinticSplineInterp;
template<typename _Ty, size_t _Options> class BisepticSplineInterp;

template<typename _Ty, size_t _Options> 
struct interpolation_traits
{ using type = BilinearInterp<_Ty, _Options>; };
template<typename _Ty> 
struct interpolation_traits<_Ty, INTERP|BICUBIC|BSPLINE>
{ using type = BicubicSplineInterp<_Ty, INTERP|BICUBIC|BSPLINE>; };
template<typename _Ty>
struct interpolation_traits<_Ty, INTERP|BIQUINTIC|BSPLINE>
{ using type = BiquinticSplineInterp<_Ty, INTERP|BIQUINTIC|BSPLINE>; };
template<typename _Ty>
struct interpolation_traits<_Ty, INTERP|BISEPTIC|BSPLINE>
{ using type = BisepticSplineInterp<_Ty, INTERP|BISEPTIC|BSPLINE>; };

template<typename _Ty, size_t _Options>
using interpolation_traits_t = typename interpolation_traits<_Ty, _Options>::type;

template<typename _Ty, typename _Derived> class InterpBase_
{
	using derived_t = _Derived;
public:
	using value_t = _Ty;
	using value_type = value_t;
	using matrix_t = types::Matrix<value_type>;

	template<typename... _Args>
	MATRICE_GLOBAL_FINL InterpBase_(const _Args&... args);

protected:
	template<typename... _Args>
	MATRICE_GLOBAL_FINL void _Bspline_coeff(const _Args& ...args);

protected:
	const value_type m_eps = value_type(1.0e-7);
	matrix_t m_coeff;
};

MATRICE_ALGS_END
#include "_base.inl"