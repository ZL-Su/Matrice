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
#include "util/_std_wrapper.h"
#include "private/_type_traits.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
struct _Blas_kernel_wrapper;
struct _Lapack_kernel_wrapper;
_DETAIL_END
using blas_kernel_t = detail::_Blas_kernel_wrapper;
using lapack_kernel_t = detail::_Lapack_kernel_wrapper;
DGE_MATRICE_END

#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
#include "..//nonfree//inl//blas_kernel_wrapper.inl"
#include "..//nonfree//inl//lapack_kernel_wrapper.inl"
#endif
