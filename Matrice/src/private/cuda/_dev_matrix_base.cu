#include "../../../include/Matrice/private/_dev_matrix_base.h"
#include "_ewise_kernels.cuh"
#include <utility>

MATRICE_DEVICE_BEGIN

template<typename _Ty>
Base_<_Ty>::pointer Base_<_Ty>::operator+(const _Myt& _other) {
	static_assert(std::is_scalar_v<_Ty>, "Template parameter _Ty must be a scalar type.");
	const int _M = *_H, _N = *_W;

	dim3 _Blocks(std::max((size_t)1, (size_t)std::ceil(_M*_N / 32.)));
	kernels::_Ewise_add<<<_Blocks, 32>>>(_Ptr, _other._Ptr, _Ptr, _M*_N);

	return (_Ptr);
}
template<typename _Ty>
Base_<_Ty>::pointer Base_<_Ty>::operator-(const _Myt& _other) {
	static_assert(std::is_scalar_v<_Ty>, "Template parameter _Ty must be a scalar type.");
	const int _M = *_H, _N = *_W;

	dim3 _Blocks(std::max((size_t)1, (size_t)std::ceil(_M*_N / 32.)));
	kernels::_Ewise_sub<<<_Blocks, 32 >>> (_Ptr, _other._Ptr, _Ptr, _M*_N);

	return (_Ptr);
}
template<typename _Ty>
Base_<_Ty>::pointer Base_<_Ty>::operator*(const _Myt& _other) {
	static_assert(std::is_scalar_v<_Ty>, "Template parameter _Ty must be a scalar type.");
	const int _M = *_H, _N = *_W;

	dim3 _Blocks(std::max((size_t)1, (size_t)std::ceil(_M*_N / 32.)));
	kernels::_Ewise_mul<<<_Blocks, 32>>> (_Ptr, _other._Ptr, _Ptr, _M*_N);

	return (_Ptr);
}
template<typename _Ty>
Base_<_Ty>::pointer Base_<_Ty>::operator/(const _Myt& _other) {
	static_assert(std::is_scalar_v<_Ty>, "Template parameter _Ty must be a scalar type.");
	const int _M = *_H, _N = *_W;

	dim3 _Blocks(std::max((size_t)1, (size_t)std::ceil(_M*_N / 32.)));
	kernels::_Ewise_div<<<_Blocks, 32>>> (_Ptr, _other._Ptr, _Ptr, _M*_N);

	return (_Ptr);
}

template class Base_<float>;
template class Base_<double>;
template class Base_<unsigned char>;
MATRICE_DEVICE_END