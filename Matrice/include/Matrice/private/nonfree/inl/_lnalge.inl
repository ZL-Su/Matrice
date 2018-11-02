#pragma once
#include "../_lnalge.h"
#ifdef __use_mkl__
#include <mkl.h>
#else
#include <fkl.h>
#endif

DGE_MATRICE_BEGIN _DETAIL_BEGIN
/**
 *\Specialization for float-type.
 */
template<> struct _Lapack_kernel_impl<float> {
	using pointer = std::add_pointer_t<float>;
	using size_type = std::tuple<int, int>;
	MATRICE_HOST_INL static int svd(pointer _A, pointer _S, pointer _Vt, const size_type& _Size) {
		auto[M, N] = _Size;
#ifdef __use_mkl__
		float* _Superb = new float[(M < N ? M : N) - 1];
		return LAPACKE_sgesvd(101, 'O', 'A',
			M, N, _A, N, _S, nullptr, 1, _Vt, N, _Superb);
#else
		return flapk::_sgesvd(_A, _S, _Vt, M, N);
#endif // __use_mkl__
	}
};

/**
 *\Specialization for double-type.
 */
template<> struct _Lapack_kernel_impl<double> {
	using pointer = std::add_pointer_t<double>;
	using size_type = std::tuple<int, int>;
	MATRICE_HOST_INL static int svd(pointer _A, pointer _S, pointer _Vt, const size_type& _Size) {
		auto[M, N] = _Size;
#ifdef __use_mkl__
		double* _Superb = new double[M < N ? M : N - 1];
		return LAPACKE_dgesvd(MKL_LAYOUT::MKL_ROW_MAJOR, 'O', 'A',
			M, N, _A, N, _S, nullptr, 1, _Vt, N, _Superb);
#else
		return flapk::_dgesvd(_A, _S, _Vt, M, N);
#endif // __use_mkl__
	}
};
_DETAIL_END DGE_MATRICE_END