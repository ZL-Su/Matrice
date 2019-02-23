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
#include "../util/_macros.h"
#include "../util/_std_wrapper.h"
//#if(defined __enable_cuda__ && !defined __disable_cuda__)
//MATRICE_PRIVATE_BEGIN
//template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
//_Ty* device_malloc(size_t& w, size_t h);
//template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
//_Ty* global_malloc(size_t size);
//template<size_t _Kind, typename _Ty>
//_Ty* device_memcpy(_Ty* hptr, _Ty* dptr, size_t w, size_t h = 1, size_t p = 1);
//template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
//void device_free(_Ty* dptr);
//template<int _Opt = 0> void _Device_sync();
//MATRICE_PRIVATE_END
//#endif