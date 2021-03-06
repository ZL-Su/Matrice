/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
**********************************************************************/
#pragma once
#include "correlation/_optim.h"

DGE_MATRICE_BEGIN
struct correlation_optimizer {
	using options = corr::detail::_Correlation_options;
	/**
	 *\brief N-th order IC-GN algorithm with bilinear interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty, uint8_t _Order = 1>
	using icgn_bilinear = corr::detail::_Corr_solver_impl<_Ty,
		bilerp_tag, corr::detail::_Alg_icgn<_Order>>;
	/**
	 *\brief N-th order IC-GN algorithm with bicubic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty, uint8_t _Order = 1>
	using icgn_bic = corr::detail::_Corr_solver_impl<_Ty,
		bicerp_tag, corr::detail::_Alg_icgn<_Order>>;
	/**
	 *\brief First order IC-GN algorithm with biquintic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty, uint8_t _Order = 1>
	using icgn_biq = corr::detail::_Corr_solver_impl<_Ty,
		biqerp_tag, corr::detail::_Alg_icgn<_Order>>;
	/**
	 *\brief First order IC-GN algorithm with biseptic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty, uint8_t _Order = 1>
	using icgn_bis = corr::detail::_Corr_solver_impl<_Ty,
		biserp_tag, corr::detail::_Alg_icgn<_Order>>;

	/**
	 *\brief First order IC-GN algorithm with bicubic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_bic_0 = corr::detail::_Corr_solver_impl<_Ty,
		bicerp_tag, corr::detail::_Alg_icgn<0>>;
	/**
	 *\brief First order IC-GN algorithm with biquintic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_biq_0 = corr::detail::_Corr_solver_impl<_Ty,
		biqerp_tag, corr::detail::_Alg_icgn<0>>;
	/**
	 *\brief First order IC-GN algorithm with biseptic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_bis_0 = corr::detail::_Corr_solver_impl<_Ty,
		biserp_tag, corr::detail::_Alg_icgn<0>>;

	/**
	 *\brief First order IC-GN algorithm with bicubic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_bic_1 = corr::detail::_Corr_solver_impl<_Ty,
		bicerp_tag, corr::detail::_Alg_icgn<1>>;
	/**
	 *\brief First order IC-GN algorithm with biquintic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_biq_1 = corr::detail::_Corr_solver_impl<_Ty,
		biqerp_tag, corr::detail::_Alg_icgn<1>>;
	/**
	 *\brief First order IC-GN algorithm with biseptic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_bis_1 = corr::detail::_Corr_solver_impl<_Ty,
		biserp_tag, corr::detail::_Alg_icgn<1>>;
};
DGE_MATRICE_END