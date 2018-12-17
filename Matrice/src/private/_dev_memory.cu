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
template<typename _Scalar, typename = typename std::enable_if<std::is_scalar<_Scalar>::value>::type>
_Scalar* device_malloc(size_t& w, size_t h)
{
	cudaError_t sts; _Scalar* dptr;
	switch (h)
	{
	case 1:
		sts = cudaMalloc(&dptr, w * sizeof(_Scalar));
		break;
	default:
		size_t pitch = 0;
		sts = cudaMallocPitch(&dptr, &pitch, w * sizeof(_Scalar), h);
		w = pitch;
		break;
	}
	if (sts != cudaSuccess) throw std::runtime_error(cudaGetErrorString(sts));
	return dptr;
}
template<typename _Scalar, typename = typename std::enable_if<std::is_scalar<_Scalar>::value>::type>
_Scalar* global_malloc(size_t N)
{
	_Scalar* dptr;
	auto sts = cudaMallocManaged(&dptr, N * sizeof(_Scalar));
	if (sts != cudaSuccess) throw std::runtime_error(cudaGetErrorString(sts));
	return dptr;
}
//<note> w is the columns, h is the rows, p is the pitch size if pitched memory is used </note>
template<typename _Scalar, int _Opt, typename = typename std::enable_if<std::is_scalar<_Scalar>::value>::type>
void device_memcpy(_Scalar* hptr, _Scalar* dptr, size_t w, size_t h = 1, size_t p = 1)
{
	cudaError_t sts; size_t hpitch = w * h * sizeof(_Scalar);
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
template<typename _Scalar, typename = typename std::enable_if<std::is_scalar<_Scalar>::value>::type>
void device_free(_Scalar* dptr) { if (dptr) cudaFree(dptr); }

#pragma region <!-- explicit intantiation -->
template int* device_malloc(size_t&, size_t);
template char* device_malloc(size_t&, size_t);
template bool* device_malloc(size_t&, size_t);
template float* device_malloc(size_t&, size_t);
template double* device_malloc(size_t&, size_t);
template unsigned char* device_malloc(size_t&, size_t);
template int* global_malloc(size_t);
template char* global_malloc(size_t);
template bool* global_malloc(size_t);
template float* global_malloc(size_t);
template double* global_malloc(size_t);
template unsigned char* global_malloc(size_t);
template void device_memcpy<int, 1>(int*, int*, size_t, size_t, size_t);
template void device_memcpy<int, 2>(int*, int*, size_t, size_t, size_t);
template void device_memcpy<char, 1>(char*, char*, size_t, size_t, size_t);
template void device_memcpy<char, 2>(char*, char*, size_t, size_t, size_t);
template void device_memcpy<bool, 1>(bool*, bool*, size_t, size_t, size_t);
template void device_memcpy<bool, 2>(bool*, bool*, size_t, size_t, size_t);
template void device_memcpy<float, 1>(float*, float*, size_t, size_t, size_t);
template void device_memcpy<float, 2>(float*, float*, size_t, size_t, size_t);
template void device_memcpy<double, 1>(double*, double*, size_t, size_t, size_t);
template void device_memcpy<double, 2>(double*, double*, size_t, size_t, size_t);
template void device_memcpy<unsigned char, 1>(unsigned char*, unsigned char*, size_t, size_t, size_t);
template void device_memcpy<unsigned char, 2>(unsigned char*, unsigned char*, size_t, size_t, size_t);
template void device_free(int*);
template void device_free(char*);
template void device_free(bool*);
template void device_free(float*);
template void device_free(double*);
template void device_free(unsigned char*);
#pragma endregion

MATRICE_PRIVATE_END

#pragma region <!-- Impl. for device_memcpy -->
template<typename T>
template<typename... _Args> MATRICE_GLOBAL
void dgelom::device::device_memcpy<T, 0>::impl(_Args ...args)
{
	return;
}
template<typename T>
template<typename... _Args> MATRICE_GLOBAL
void dgelom::device::device_memcpy<T, 1>::impl(_Args ...args)
{
	dgelom::privt::device_memcpy<T, option>(args...);
}
template void dgelom::device::device_memcpy<uchar, 1>
::impl(pointer, pointer, size_t, size_t, size_t);
template void dgelom::device::device_memcpy<float, 1>
::impl(pointer, pointer, size_t, size_t, size_t);
template void dgelom::device::device_memcpy<double, 1>
::impl(pointer, pointer, size_t, size_t, size_t);
template<typename T>
template<typename... _Args> MATRICE_GLOBAL
void dgelom::device::device_memcpy<T, 2>::impl(_Args ...args)
{
	dgelom::privt::device_memcpy<T, option>(args...);
}
template void dgelom::device::device_memcpy<uchar, 2>
::impl(pointer, pointer, size_t, size_t, size_t);
template void dgelom::device::device_memcpy<float, 2>
::impl(pointer, pointer, size_t, size_t, size_t);
template void dgelom::device::device_memcpy<double, 2>
::impl(pointer, pointer, size_t, size_t, size_t);
#pragma endregion

#endif