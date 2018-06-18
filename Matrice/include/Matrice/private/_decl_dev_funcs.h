/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
**************************************************************************/
#pragma once
#include <type_traits>
#include "../util/_macros.h"
#if(defined __enable_cuda__ && !defined __disable_cuda__)
MATRICE_PRIVATE_BEGIN
template<typename _Scalar, typename = typename std::enable_if<std::is_scalar<_Scalar>::value>::type>
_Scalar* device_malloc(std::size_t& w, std::size_t h);
template<typename _Scalar, typename = typename std::enable_if<std::is_scalar<_Scalar>::value>::type>
_Scalar* global_malloc(std::size_t size);
template<typename _Scalar, int _Opt, typename = typename std::enable_if<std::is_scalar<_Scalar>::value>::type>
_Scalar* device_memcpy(_Scalar* hptr, _Scalar* dptr, size_t w, size_t h = 1, size_t p = 1);
template<typename _Scalar, typename = typename std::enable_if<std::is_scalar<_Scalar>::value>::type>
void device_free(_Scalar* dptr);
template<int _Opt = 0> void _Device_sync();
MATRICE_PRIVATE_END
#endif