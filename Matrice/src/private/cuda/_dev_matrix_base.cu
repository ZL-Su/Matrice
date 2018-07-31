#include "../../../include/Matrice/private/_dev_matrix_base.h"
#include "_ewise_kernels.cuh"

MATRICE_DEVICE_BEGIN

template<typename _Ty, typename _Derived>
_Derived& Base_<_Ty, _Derived>::operator+(const _Derived& _other) {
	const auto _M = *_other._H, _N = *_other._W;

	dim3 _Blocks(static_cast<size_t>(ceil(_M*_N / 512.)));
	kernels::_Ewise_plus<<<_Blocks, 512>>>(_Ptr, _other._Ptr, _Ptr, _M*_N);

	return (*static_cast<_Derived*>(this));
}
template<typename _Ty> using DMatrix = types::Matrix_<_Ty, -1, -1>;
DMatrix<float>& Base_<float, DMatrix<float>>::operator+(const DMatrix<float>&);
DMatrix<double>& Base_<double, DMatrix<double>>::operator+(const DMatrix<double>&);

MATRICE_DEVICE_END