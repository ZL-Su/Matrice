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
#include "../private/_memory.h"
#if(defined MATRICE_ENABLE_CUDA && !defined __disable_cuda__)
#include <cuda_runtime.h>
MATRICE_DEVICE_BEGIN
enum { LINEAR = 8000, PITCHED = 8001, ARRTARR = 8002, FROMARR = 8003, TOARRAY = 8004, };
///<!-- generic first-type base for device operators -->
template<typename T> struct device_base_v1
{
	using value_t = T; 
	using pointer = value_t*;
	using hptr_t = pointer;
	using dptr_t = pointer;
	using const_hptr = const hptr_t;
	using const_dptr = const dptr_t;
};
///<!-- generic first-type class for device memory operation-->
template<typename T, size_t _CpyKind>
struct device_memcpy { enum { option = _CpyKind}; };
template<typename T> struct device_memcpy<T, 0> : device_base_v1<T>
{
	using typename device_base_v1<T>::pointer;
	enum {option = 0};
	template<typename... _Args>
	static MATRICE_GLOBAL void impl(_Args... args);
};
template<typename T> struct device_memcpy<T, 1> : device_base_v1<T>
{
	using typename device_base_v1<T>::pointer;
	enum { option = 1 };
	template<typename... _Args>
	static MATRICE_GLOBAL void impl(_Args... args) {
		privt::device_memcpy<option>(args...);
	}
};
template<typename T> struct device_memcpy<T, 2> : device_base_v1<T>
{
	using typename device_base_v1<T>::pointer;
	enum { option = 2 };
	template<typename... _Args>
	static MATRICE_GLOBAL void impl(_Args... args) {
		privt::device_memcpy<option>(args...);
	}
};

MATRICE_DEVICE_END
#endif

