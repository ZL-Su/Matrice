/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2024, Zhilong(Dgelom) Su, all rights reserved.

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
#include "Matrice/private/_memory.h"

#ifdef MATRICE_ENABLE_CUDA
#include <cuda_runtime.h>
#pragma warning(disable: 4715 4661 4224 4267 4244 4819 4199)

MATRICE_PRIVATE_BEGIN
//<note> w is the columns, h is the rows </note>
template<typename _Ty, typename>
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
template<typename _Ty, typename>
_Ty* global_malloc(size_t N) {
	_Ty* dptr; auto sts = cudaMallocManaged(&dptr, N * sizeof(_Ty));
	if (sts != cudaSuccess)
		throw std::runtime_error(cudaGetErrorString(sts));
	return dptr;
}
//<note> w is the columns, h is the rows, p is the pitch size if pitched memory is used </note>
template<size_t _Kind, typename _Ty, typename>
void device_memcpy(_Ty* hptr, _Ty* dptr, size_t w, size_t h, size_t p) {
	cudaError_t sts; size_t hpitch = w * h * sizeof(_Ty);
	if (_Kind == ::cudaMemcpyHostToDevice) {
		if (p == 1) {
			sts = cudaMemcpy(dptr, hptr, hpitch, cudaMemcpyHostToDevice);
		}
		else {
			sts = cudaMemcpy2D(dptr, p, hptr, hpitch, hpitch, 1, cudaMemcpyHostToDevice);
		}
	}
	if (_Kind == ::cudaMemcpyDeviceToHost) {
		if (p == 1) {
			sts = cudaMemcpy(hptr, dptr, hpitch, ::cudaMemcpyDeviceToHost);
		}
		else {
			sts = cudaMemcpy2D(hptr, hpitch, dptr, p, hpitch, 1, ::cudaMemcpyDeviceToHost);
		}
	}
	if (sts != cudaSuccess) throw std::runtime_error(cudaGetErrorString(sts));
}
template<typename _Ty, typename>
void device_free(_Ty* dptr) { 
	if (dptr) {
		auto sts = cudaFree(dptr);
		if(sts != cudaSuccess)
			throw std::runtime_error(cudaGetErrorString(sts));
	}
}

#pragma region <!-- explicit intantiation -->
#define _DEVMALLOC(T) \
template T* device_malloc(size_t&, size_t); \
template T* global_malloc(size_t); \
template void device_free(T*);

#define _DEVMEMCPY(T) \
template void device_memcpy<1,T,void>(T*, T*, size_t, size_t, size_t); \
template void device_memcpy<2,T,void>(T*, T*, size_t, size_t, size_t);

_DEVMALLOC(int)            _DEVMALLOC(char)
_DEVMALLOC(unsigned char)  _DEVMALLOC(bool)
_DEVMALLOC(float)          _DEVMALLOC(double)
_DEVMALLOC(size_t)
_DEVMEMCPY(int)            _DEVMEMCPY(char)
_DEVMEMCPY(bool)           _DEVMEMCPY(unsigned char)
_DEVMEMCPY(float)          _DEVMEMCPY(double)
_DEVMEMCPY(size_t)

#pragma endregion

MATRICE_PRIVATE_END

#undef _DEVMALLOC
#undef _DEVMEMCPY
#endif