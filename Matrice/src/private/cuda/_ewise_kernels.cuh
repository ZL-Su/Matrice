#pragma once

#include "../../../include/Matrice/util/_macros.h"

#ifdef __enable_cuda__
MATRICE_DEVICE_BEGIN namespace kernels {
template<typename _Scalar> __global__ void
_Ewise_reduce(const _Scalar* _Data, size_t _Size, _Scalar& _Ret) {
	_Ret = _Scalar(0);
	
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;
	if (_Idx < _Size) _Ret += _Data[_Idx];
}

template<typename _Scalar> __global__ void
_Ewise_plus(const _Scalar* _X, const _Scalar* _Y, _Scalar* _Z, size_t _Size) {
	const auto _Idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (_Idx < _Size) _Z[_Idx] = _X[_Idx] + _Y[_Idx];
}
}MATRICE_DEVICE_END
#endif