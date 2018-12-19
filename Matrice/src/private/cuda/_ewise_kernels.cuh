#pragma once
#include <type_traits>
#include "../../../include/Matrice/util/_macros.h"

#ifdef __enable_cuda__
#include <cublas.h>
#include <thrust\device_vector.h>

MATRICE_DEVICE_BEGIN namespace kernels {
template<typename _Scalar, typename = std::enable_if_t<std::is_scalar_v<_Scalar>>> 
__global__ void _Ewise_reduce(const _Scalar* _Data, std::size_t _Size, _Scalar& _Ret) {
	_Ret = _Scalar(0);
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) _Ret += _Data[_Idx];
}

/**
 * $Z_i = X_i + Y_i$
 */
template<typename _Scalar> __global__ void
_Ewise_add(const _Scalar* _X, const _Scalar* _Y, _Scalar* _Z, std::size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] + _Y[_Idx];
	}
}
/**
 * $Z_i = X_i - Y_i$
 */
template<typename _Scalar> __global__ void
_Ewise_sub(const _Scalar* _X, const _Scalar* _Y, _Scalar* _Z, std::size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] - _Y[_Idx];
	}
}
/**
 * $Z_i = X_i * Y_i$
 */
template<typename _Scalar> __global__ void
_Ewise_mul(const _Scalar* _X, const _Scalar* _Y, _Scalar* _Z, std::size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] * _Y[_Idx];
	}
}
/**
 * $Z_i = X_i / Y_i$
 */
template<typename _Scalar> __global__ void
_Ewise_div(const _Scalar* _X, const _Scalar* _Y, _Scalar* _Z, std::size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) {
		_Z[_Idx] = _X[_Idx] / _Y[_Idx];
	}
}
}MATRICE_DEVICE_END
#endif