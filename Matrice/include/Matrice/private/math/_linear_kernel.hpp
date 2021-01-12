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
#include "_primitive_funcs.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ptr> MATRICE_GLOBAL
/// <summary>
/// Cholesky decomposition.
/// </summary>
/// <param name="data"></param>
/// <param name="n"></param>
/// <returns></returns>
int _Linear_spd_kernel(_Ptr data, size_t n) noexcept;

template<typename _Ptr> MATRICE_GLOBAL
/// <summary>
/// Matrix inversion with Cholesky decomposition.
/// </summary>
/// <param name="data"></param>
/// <param name="inv"></param>
/// <param name="n"></param>
/// <returns></returns>
void _Linear_ispd_kernel(_Ptr data, _Ptr inv, size_t n) noexcept;

template<typename _Ptr> MATRICE_GLOBAL
/// <summary>
/// Solve Cholesky factorized linear system 'L^T*L*x = x'.
/// </summary>
/// <param name="n"></param>
/// <param name="l"></param>
/// <param name="x"></param>
/// <param name="stride"></param>
/// <returns></returns>
void _Linear_spd_bwd(size_t n, _Ptr l, _Ptr x, int stride=1) noexcept;

template<typename _Ptr> MATRICE_GLOBAL
/// <summary>
/// Perform LU factorization.
/// </summary>
/// <param name="n"></param>
/// <param name="data"></param>
/// <param name="idx"></param>
/// <returns></returns>
int _Linear_lud_kernel(size_t n, _Ptr data, int* idx) noexcept;

template<typename _Ptr> MATRICE_GLOBAL
/// <summary>
/// Solve LU factorized linear system 'lu*x = x'.
/// </summary>
/// <param name="n"></param>
/// <param name="lu"></param>
/// <param name="x"></param>
/// <param name="stride"></param>
/// <returns></returns>
void _Linear_lud_sv(size_t n, _Ptr lu, _Ptr x, int stride = 1) noexcept;

_DETAIL_END
DGE_MATRICE_END