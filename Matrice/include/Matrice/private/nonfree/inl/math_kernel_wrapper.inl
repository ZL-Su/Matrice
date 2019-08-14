/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include "../../math/kernel_wrapper.hpp"
#include "../../math/_config.h"

DGE_MATRICE_BEGIN
_INTERNAL_BEGIN
template<typename _It>
decltype(auto) _Blas_asum(const _It x, size_t n, int inc) {};
template<typename _Ty>
decltype(auto) _Blas_axpy(const _Ty a, const _Ty* x, _Ty* y, size_t n, int incx, int incy) {};
_INTERNAL_END

_DETAIL_BEGIN
struct _Blas_kernel_wrapper {
	template<class _Aty, typename _Ret = typename _Aty::value_type>
	static MATRICE_HOST_INL _Ret reduce(const _Aty& a, int inc = 1) noexcept {
		return internal::_Blas_asum(a.data(), a.size(), inc);
	}

	template<class _Xty, class _Yty=_Xty, 
		typename = enable_if_t<is_same_v<typename _Xty::value_type, typename _Yty::value_type>>>
	static MATRICE_HOST_INL _Yty& axpy(const _Xty& x, _Yty& y, typename _Xty::value_type a=1) noexcept {
		internal::_Blas_axpy(a, x.data(), y.data(), min(x.size(), y.size()), 1, 1);
		return (y);
	}
};
_DETAIL_END

_INTERNAL_BEGIN
template<>
decltype(auto) _Blas_asum(const float* x, size_t n, int inc) {
	return cblas_sasum(MKL_INT(n), x, MKL_INT(inc));
}
template<>
decltype(auto) _Blas_asum(const double* x, size_t n, int inc) {
	return cblas_dasum(MKL_INT(n), x, MKL_INT(inc));
}

template<>
decltype(auto) _Blas_axpy(const float a, const float* x, float* y, size_t n, int incx, int incy) {
	return cblas_saxpy(n, a, x, incx, y, incy);
}
template<>
decltype(auto) _Blas_axpy(const double a, const double* x, double* y, size_t n, int incx, int incy) {
	return cblas_daxpy(n, a, x, incx, y, incy);
}
_INTERNAL_END
DGE_MATRICE_END
