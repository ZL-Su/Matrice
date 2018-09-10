#include "../../../include/Matrice/private/_dev_matrix_base.h"
#include "_ewise_kernels.cuh"

MATRICE_DEVICE_BEGIN

template<typename arith_type>
Base_<arith_type>::pointer Base_<arith_type>::operator+(const pointer _other) {
	const int _M = *_H, _N = *_W;

	dim3 _Blocks(ceil<size_t>(_M*_N / 32.));
	kernels::_Ewise_plus<<<_Blocks, 32>>>(_Ptr, _other, _Ptr, _M*_N);

	return (_Ptr);
}

template class Base_<float>;
template class Base_<double>;
template class Base_<unsigned char>;
MATRICE_DEVICE_END