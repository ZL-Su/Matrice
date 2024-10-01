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
#include <utility>
#include "Matrice/private/_dev_matrix_base.h"

#ifdef MATRICE_ENABLE_CUDA
#include "_ewise_kernels.cuh"

MATRICE_DEVICE_BEGIN

#define _KERNEL_CONFIG <<_Blocks, threads_per_block>>

template<typename _Ty>
Base_<_Ty>::pointer Base_<_Ty>::operator+(const _Myt& _other) {
	static_assert(std::is_scalar_v<_Ty>, "Template parameter _Ty must be a scalar type.");
	const auto N = *_H**_W;

	dim3 _Blocks((N + threads_per_block - 1) / threads_per_block);
	kernels::_Ewise_add<_KERNEL_CONFIG>(_Ptr, _other._Ptr, _Ptr, N);

	return (_Ptr);
}
template<typename _Ty>
Base_<_Ty>::pointer Base_<_Ty>::operator-(const _Myt& _other) {
	static_assert(std::is_scalar_v<_Ty>, "Template parameter _Ty must be a scalar type.");
	const auto N = *_H**_W;

	dim3 _Blocks((N + threads_per_block - 1) / threads_per_block);
	kernels::_Ewise_sub<_KERNEL_CONFIG>(_Ptr, _other._Ptr, _Ptr, N);

	return (_Ptr);
}
template<typename _Ty>
Base_<_Ty>::pointer Base_<_Ty>::operator*(const _Myt& _other) {
	static_assert(std::is_scalar_v<_Ty>, "Template parameter _Ty must be a scalar type.");
	const auto N = *_H**_W;

	dim3 _Blocks((N + threads_per_block - 1) / threads_per_block);
	kernels::_Ewise_mul<_KERNEL_CONFIG>(_Ptr, _other._Ptr, _Ptr, N);

	return (_Ptr);
}
template<typename _Ty>
Base_<_Ty>::pointer Base_<_Ty>::operator/(const _Myt& _other) {
	static_assert(std::is_scalar_v<_Ty>, "Template parameter _Ty must be a scalar type.");
	const auto N = *_H**_W;

	dim3 _Blocks((N+threads_per_block-1)/ threads_per_block);
	kernels::_Ewise_div<_KERNEL_CONFIG>(_Ptr, _other._Ptr, _Ptr, N);

	return (_Ptr);
}

template<typename _Ty>
Base_<_Ty>::value_t Base_<_Ty>::reduce() {
	auto _Ret = value_t(0);

	const auto N = *_H**_W;
	dim3 _Blocks((N + threads_per_block - 1) / threads_per_block);
	kernels::_Ewise_reduce<_KERNEL_CONFIG>(_Ret, _Ptr, N);
	_Sync_impl();

	return (_Ret);
}

template<typename _Ty>
typename Base_<_Ty>::_Mydt operator+(const typename Base_<_Ty>::value_t _Left, const typename Base_<_Ty>::_Mydt& _Right) {
	typename Base_<_Ty>::_Mydt _Ret(_Right.rows(), _Right.cols());

	const auto _N = _Ret.size();
	dim3 _Blocks((_N + threads_per_block - 1) / threads_per_block);
	kernels::_Ewise_add<_KERNEL_CONFIG>(_Right.data(), _Left, _Ret.data(), _N);

	return (_Ret);
}

///<brief> explicit instantiations </brief>
template class Base_<float>;
template class Base_<double>;
template class Base_<unsigned char>;

#undef _KERNEL_CONFIG
MATRICE_DEVICE_END
#endif