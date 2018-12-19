#include "../../../include/Matrice/private/_dev_matrix_base.h"
#include "_ewise_kernels.cuh"
#include <utility>

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
Base_<_Ty>::value_t Base_<_Ty>::reduce() const {
	auto _Ret = value_t(0);

	const auto N = *_H**_W;
	dim3 _Blocks((N + threads_per_block - 1) / threads_per_block);
	kernels::_Ewise_reduce<_KERNEL_CONFIG>(_Ret, _Ptr, N);
	_Sync_impl();

	return (_Ret);
}


///<brief> explicit instantiations </brief>
template class Base_<float>;
template class Base_<double>;
template class Base_<unsigned char>;

#undef _KERNEL_CONFIG
MATRICE_DEVICE_END