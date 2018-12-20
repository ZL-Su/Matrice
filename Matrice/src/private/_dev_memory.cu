/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#include <complex>
#include <stdexcept>
#include "../../include/Matrice/util/_macros.h"
#include "../../include/Matrice/private/_devops.h"
#include "../../include/Matrice/private/_unified_memory.h"

#if (defined __enable_cuda__ && !defined __disable_cuda__)
#include <cuda_runtime.h>
#pragma warning(disable: 4715 4661 4224 4267 4244 4819 4199)

using std::size_t;
using std::complex;
using uchar = unsigned char;
MATRICE_PRIVATE_BEGIN
//<note> w is the columns, h is the rows </note>
template<typename _Ty, typename = std::enable_if_t<std::is_scalar_v<_Ty>>>
_Ty* device_malloc(size_t& w, size_t h) {
	cudaError_t sts; _Ty* dptr;
	switch (h)
	{
	case 1:
		sts = cudaMalloc(&dptr, w * sizeof(_Ty));
		break;
	default:
		size_t pitch = 0;
		sts = cudaMallocPitch(&dptr, &pitch, w * sizeof(_Ty), h);
		w = pitch;
		break;
	}
	if (sts != cudaSuccess) throw std::runtime_error(cudaGetErrorString(sts));
	return dptr;
}
template<typename _Ty, typename = std::enable_if_t<std::is_scalar_v<_Ty>>>
_Ty* global_malloc(size_t N) {
	_Ty* dptr; auto sts = cudaMallocManaged(&dptr, N * sizeof(_Ty));
	if (sts != cudaSuccess)
		throw std::runtime_error(cudaGetErrorString(sts));
	return dptr;
}
//<note> w is the columns, h is the rows, p is the pitch size if pitched memory is used </note>
template<typename _Ty, int _Opt, typename = std::enable_if_t<std::is_scalar_v<_Ty>>>
void device_memcpy(_Ty* hptr, _Ty* dptr, size_t w, size_t h = 1, size_t p = 1) {
	cudaError_t sts; size_t hpitch = w * h * sizeof(_Ty);
	if (_Opt == ::cudaMemcpyHostToDevice) {
		if (p == 1) {
			sts = cudaMemcpy(dptr, hptr, hpitch, cudaMemcpyHostToDevice);
		}
		else {
			sts = cudaMemcpy2D(dptr, p, hptr, hpitch, hpitch, 1, cudaMemcpyHostToDevice);
		}
	}
	if (_Opt == ::cudaMemcpyDeviceToHost) {
		if (p == 1) {
			sts = cudaMemcpy(hptr, dptr, hpitch, ::cudaMemcpyDeviceToHost);
		}
		else {
			sts = cudaMemcpy2D(hptr, hpitch, dptr, p, hpitch, 1, ::cudaMemcpyDeviceToHost);
		}
	}
	if (sts != cudaSuccess) throw std::runtime_error(cudaGetErrorString(sts));
}
template<typename _Ty, typename = std::enable_if_t<std::is_scalar_v<_Ty>>>
void device_free(_Ty* dptr) { 
	if (dptr) {
		auto sts = cudaFree(dptr);
		if(sts != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(sts));
	}
}

#pragma region <!-- explicit intantiation -->
#define _DEVMALLOC(type) \
template type* device_malloc(size_t&, size_t); \
template type* global_malloc(size_t); \
template void device_free(type*);
#define _DEVMEMCPY(type, op) \
template void device_memcpy<type, op>(type*, type*, size_t, size_t, size_t);

_DEVMALLOC(int)        _DEVMALLOC(char)
_DEVMALLOC(uchar)      _DEVMALLOC(bool)
_DEVMALLOC(float)      _DEVMALLOC(double)
_DEVMEMCPY(int, 1)     _DEVMEMCPY(int, 2)
_DEVMEMCPY(char, 1)    _DEVMEMCPY(char, 2)
_DEVMEMCPY(bool, 1)    _DEVMEMCPY(bool, 2)
_DEVMEMCPY(uchar, 1)   _DEVMEMCPY(uchar, 2)
_DEVMEMCPY(float, 1)   _DEVMEMCPY(float, 2)
_DEVMEMCPY(double, 1)  _DEVMEMCPY(double, 2)

#pragma endregion

MATRICE_PRIVATE_END

#pragma region <!-- Impl. for device_memcpy -->
#define _DEVMEMCPY_IMPL(type, option) \
template void dgelom::device::device_memcpy<type, option> \
			::impl(pointer, pointer, size_t, size_t, size_t);

template<typename T>
template<typename... _Args> MATRICE_GLOBAL
void dgelom::device::device_memcpy<T, 0>::impl(_Args ...args) {
	return;
}
template<typename T>
template<typename... _Args> MATRICE_GLOBAL
void dgelom::device::device_memcpy<T, 1>::impl(_Args ...args) {
	dgelom::privt::device_memcpy<T, option>(args...);
}
template<typename T>
template<typename... _Args> MATRICE_GLOBAL
void dgelom::device::device_memcpy<T, 2>::impl(_Args ...args) {
	dgelom::privt::device_memcpy<T, option>(args...);
}

_DEVMEMCPY_IMPL(uchar, 1)
_DEVMEMCPY_IMPL(float, 1)
_DEVMEMCPY_IMPL(double, 1)
_DEVMEMCPY_IMPL(uchar,  2)
_DEVMEMCPY_IMPL(float,  2)
_DEVMEMCPY_IMPL(double, 2)
#pragma endregion

#undef _DEVMALLOC
#undef _DEVMEMCPY
#undef _DEVMEMCPY_IMPL
#endif