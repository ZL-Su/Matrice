#include "../../../include/Matrice/private/_dev_matrix_base.h"
#include "_ewise_kernels.cuh"

MATRICE_DEVICE_BEGIN
template<typename _Ty, typename _Derived>
_Derived& Base_<_Ty, _Derived>::operator+(const Base_& _other) {
	const auto _M = *_other._H, _N = *_other._W;

	dim3 _Blocks(static_cast<size_t>(ceil(_M*_N / 512.)));
	kernels::_Ewise_plus<<<_Blocks, 512>>>(_Ptr, _other._Ptr, _Ptr, _M*_N);

	return (*static_cast<_Derived*>(this));
}

types::Matrix_<float, -1, -1>& Base_<float, types::Matrix_<float, -1, -1>>::operator+(const Base_<value_t, types::Matrix_<value_t, -1, -1>>&);
types::Matrix_<double, -1, -1>& Base_<double, types::Matrix_<double, -1, -1>>::operator+(const Base_<value_t, types::Matrix_<value_t, -1, -1>>&);
MATRICE_DEVICE_END