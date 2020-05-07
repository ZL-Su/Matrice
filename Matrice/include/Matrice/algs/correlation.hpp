/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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
#include "correlation/_optim.h"

DGE_MATRICE_BEGIN
struct correlation_optimizer {
	using options = algs::detail::corr::_Correlation_options;

	/**
	 *\brief N-th order IC-GN alg. with bicubic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty, uint8_t _Order = 1>
	using icgn_bic = algs::detail::corr::_Corr_solver_impl<_Ty,
		bicerp_tag, algs::detail::corr::_Alg_icgn<_Order>>;
	/**
	 *\brief 1th order IC-GN alg. with biquintic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty, uint8_t _Order = 1>
	using icgn_biq = algs::detail::corr::_Corr_solver_impl<_Ty,
		biqerp_tag, algs::detail::corr::_Alg_icgn<_Order>>;
	/**
	 *\brief 1th order IC-GN alg. with biseptic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty, uint8_t _Order = 1>
	using icgn_bis = algs::detail::corr::_Corr_solver_impl<_Ty,
		biserp_tag, algs::detail::corr::_Alg_icgn<_Order>>;

	/**
	 *\brief 1th order IC-GN alg. with bicubic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_bic_0 = algs::detail::corr::_Corr_solver_impl<_Ty,
		bicerp_tag, algs::detail::corr::_Alg_icgn<0>>;
	/**
	 *\brief 1th order IC-GN alg. with biquintic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_biq_0 = algs::detail::corr::_Corr_solver_impl<_Ty,
		biqerp_tag, algs::detail::corr::_Alg_icgn<0>>;
	/**
	 *\brief 1th order IC-GN alg. with biseptic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_bis_0 = algs::detail::corr::_Corr_solver_impl<_Ty,
		biserp_tag, algs::detail::corr::_Alg_icgn<0>>;

	/**
	 *\brief 1th order IC-GN alg. with bicubic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_bic_1 = algs::detail::corr::_Corr_solver_impl<_Ty,
		bicerp_tag, algs::detail::corr::_Alg_icgn<1>>;
	/**
	 *\brief 1th order IC-GN alg. with biquintic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_biq_1 = algs::detail::corr::_Corr_solver_impl<_Ty,
		biqerp_tag, algs::detail::corr::_Alg_icgn<1>>;
	/**
	 *\brief 1th order IC-GN alg. with biseptic spline interpolation.
	 *\param <_Ty> must be a scalar type of float or double.
	 */
	template<typename _Ty>
	using icgn_bis_1 = algs::detail::corr::_Corr_solver_impl<_Ty,
		biserp_tag, algs::detail::corr::_Alg_icgn<1>>;
};
DGE_MATRICE_END