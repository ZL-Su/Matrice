/*  *************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include <tuple>
#include "../_type_traits.h"

DGE_MATRICE_BEGIN _DETAIL_BEGIN

/**
 *\brief blas kernels for matrix and vector linear operation
 *\param <_Ty> should be a scalar type.
 */
template<typename _Ty, MATRICE_ENABLE_IF(is_floating_point_v<_Ty>)> 
struct _Blas_kernel_impl {
	static_assert("Oops, unsupported data type _Ty in _Blas_kernel_impl<_Ty, void>.");

	template<typename... _Args> static constexpr auto dot(const _Args&...) {}
	template<typename... _Args> static constexpr auto mul(const _Args&...) {}
	template<typename... _Args> static constexpr auto gemm(const _Args&...) {}
};

/**
 *\brief lapack kernels for matrix and vector linear operation
 *\param <_Ty> should be a scalar type.
 */
template<typename _Ty, MATRICE_ENABLE_IF(is_floating_point_v<_Ty>)>
struct _Lapack_kernel_impl {
	static_assert("Oops, unsupported data type _Ty in _Lapack_kernel_impl<_Ty, void>.");

	template<typename... _Args> static constexpr auto svd(const _Args&...) {}
	template<typename... _Args> static constexpr auto spd(const _Args&...) {}
	template<typename... _Args> static constexpr auto lud(const _Args&...) {}
	template<typename... _Args> static constexpr auto slv(const _Args&...) {}
};

/**
 *\brief lapack backward for solving linear eqautions
 *\param <_Tag> a linear solver type in solver_type::AUTO/LUF/CHD/QRD/SVD.
 */
template<solver_type _Tag> struct _Lapack_backward_impl {
	template<typename... _Args> 
	static constexpr auto eval(const _Args&...) {}
};

_DETAIL_END
template<typename _Ty>
using blas_kernel = detail::_Blas_kernel_impl<_Ty>;
template<typename _Ty>
using lapack_kernel = detail::_Lapack_kernel_impl<_Ty>;
template<solver_type _Tag>
using lapack_backward = detail::_Lapack_backward_impl<_Tag>;
DGE_MATRICE_END

#include "inl\_lnalge.inl"