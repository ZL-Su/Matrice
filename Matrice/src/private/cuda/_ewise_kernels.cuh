/**************************************************************************
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
**************************************************************************/
#pragma once
#include <type_traits>
#include "Matrice/util/_macros.h"

#ifdef MATRICE_ENABLE_CUDA
#include <cublas.h>
#include <thrust\device_vector.h>

MATRICE_DEVICE_BEGIN namespace kernels {

/**
 * summation over _Data[0 : _Size - 1]
 */
template<typename _Ty> __global__ 
void _Ewise_reduce(_Ty& _Ret, const _Ty* _Data, size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) _Ret = _Ret + _Data[_Idx];
}

/**
 * $Z_i = X_i + Y_i$
 */
template<typename _Ty> __global__ 
void _Ewise_add(const _Ty* _X, const _Ty* _Y, _Ty* _Z, size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] + _Y[_Idx];
	}
}
template<typename _Ty> __global__
void _Ewise_add(const _Ty* _X, const _Ty _Y, _Ty* _Z, size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] + _Y;
	}
}
/**
 * $Z_i = X_i - Y_i$
 */
template<typename _Ty> __global__ 
void _Ewise_sub(const _Ty* _X, const _Ty* _Y, _Ty* _Z, size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] - _Y[_Idx];
	}
}
/**
 * $Z_i = X_i * Y_i$
 */
template<typename _Ty> __global__ 
void _Ewise_mul(const _Ty* _X, const _Ty* _Y, _Ty* _Z, size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] * _Y[_Idx];
	}
}
template<typename _Ty> __global__
void _Ewise_mul(const _Ty* _X, const _Ty _Y, _Ty* _Z, size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] * _Y;
	}
}
/**
 * $Z_i = X_i / Y_i$
 */
template<typename _Ty> __global__ 
void _Ewise_div(const _Ty* _X, const _Ty* _Y, _Ty* _Z, size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] / _Y[_Idx];
	}
}
}MATRICE_DEVICE_END
#endif