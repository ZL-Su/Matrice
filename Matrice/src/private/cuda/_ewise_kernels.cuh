#pragma once
#include <type_traits>
#include "../../../include/Matrice/util/_macros.h"

#ifdef __enable_cuda__
#include <cublas.h>
#include <thrust\device_vector.h>

MATRICE_DEVICE_BEGIN namespace kernels {

/**
 * summation over _Data[0 : _Size - 1]
 */
template<typename _Ty> __global__ 
void _Ewise_reduce(_Ty& _Ret, const _Ty* _Data, std::size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) _Ret = _Ret + _Data[_Idx];
}

/**
 * $Z_i = X_i + Y_i$
 */
template<typename _Ty> __global__ 
void _Ewise_add(const _Ty* _X, const _Ty* _Y, _Ty* _Z, std::size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] + _Y[_Idx];
	}
}
/**
 * $Z_i = X_i - Y_i$
 */
template<typename _Ty> __global__ 
void _Ewise_sub(const _Ty* _X, const _Ty* _Y, _Ty* _Z, std::size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] - _Y[_Idx];
	}
}
/**
 * $Z_i = X_i * Y_i$
 */
template<typename _Ty> __global__ 
void _Ewise_mul(const _Ty* _X, const _Ty* _Y, _Ty* _Z, std::size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] * _Y[_Idx];
	}
}
/**
 * $Z_i = X_i / Y_i$
 */
template<typename _Ty> __global__ 
void _Ewise_div(const _Ty* _X, const _Ty* _Y, _Ty* _Z, std::size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] / _Y[_Idx];
	}
}
}MATRICE_DEVICE_END
#endif