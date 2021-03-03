/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for
more detail.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once
#include "util/_macros.h"
#include "private/_type_traits.h"

DGE_MATRICE_BEGIN
// Tag base
struct interpolation_tag {};

// Tag, \sa [l]inear int[erp]olation
struct lerp_tag   :interpolation_tag {};

// Tag, \sa [bil]inear int[erp]olation
struct bilerp_tag :interpolation_tag {};

// Tag, \sa [cub]ic [con]volution int[erp]olation
struct cubconerp_tag :interpolation_tag {};

// Tag, \sa [bic]ubic int[erp]olation
struct bicerp_tag :interpolation_tag {};

// Tag, \sa [biq]intic int[erp]olation
struct biqerp_tag :interpolation_tag {};

// Tag, \sa [bis]eptic int[erp]olation
struct biserp_tag :interpolation_tag {};

// Tag, \sa 2d [m]ultilevel bicubic interpolation
struct mbicerp_2d_tag :interpolation_tag {
	static constexpr size_t dimension = 2;
	static constexpr size_t max_levels = 8;
};

// Tag, \sa 3d [m]ultilevel bicubic interpolation
struct mbicerp_3d_tag :interpolation_tag {
	static constexpr size_t dimension = 3;
	static constexpr size_t max_levels = 8;
};

// Forward declarations
namespace algs {
template<typename _Ty, typename _Tag> class _Interpolation_wrapper;
template<typename _Ty, typename _Tag> class auto_interp_dispatcher;
template<typename _Ty, unsigned _N>   class multilevel_bicerp_approx;
}

/// <summary>
/// ALIAS class template for b-spline interpolator.
/// </summary>
/// <typeparam name="_Ty">data type</typeparam>
/// <typeparam name="_Tag">interpolation algorithm tag</typeparam>
template<
	typename _Ty, 
	typename _Tag = bicerp_tag, 
	MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
using interpolation = algs::_Interpolation_wrapper<_Ty, _Tag>;

/// <summary>
/// ALIAS class template for interpolator dispatcher.
/// </summary>
/// <typeparam name="_Ty">data type</typeparam>
/// <typeparam name="_Tag">interpolation algorithm tag</typeparam>
template<
	typename _Ty,
	typename _Tag = bicerp_tag,
	MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
using interp_alg_adapter = algs::auto_interp_dispatcher<_Ty, _Tag>;

/**
 *\brief 1-dimensional scattered data interpolation.
 *\param <_Ty> data type, must be float or double.
 */
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
using multilevel_bicerp_1d = algs::multilevel_bicerp_approx<_Ty, 1>;

/**
 *\brief 2-dimensional scattered data interpolation.
 *\param <_Ty> data type, must be float or double.
 */
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
using multilevel_bicerp_2d = algs::multilevel_bicerp_approx<_Ty, 2>;

/**
 *\brief 3-dimensional scattered data interpolation.
 *\param <_Ty> data type, must be float or double.
 */
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
using multilevel_bicerp_3d = algs::multilevel_bicerp_approx<_Ty, 3>;
DGE_MATRICE_END