#pragma once
#include "../_lnalge.h"
#ifdef __use_mkl__
#include <mkl.h>
#else
#include <fkl.h>
#endif

DGE_MATRICE_BEGIN _DETAIL_BEGIN
template<typename _Ty> struct _Lapack_kernel_impl_base {
	using pointer = std::add_pointer_t<_Ty>;
	using size_type = std::tuple<int, int>;
	using plview_type = std::tuple<int, int, pointer>;
};
/**
 *\Specialization for float-type.
 */
template<> struct _Lapack_kernel_impl<float> : _Lapack_kernel_impl_base<float> {
	/**
	 * \computes singular value decomposition
	 * \Output: $_A := U, _S := \Sigma, _Vt := V^T$
	 */
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

	/**
	 * \computes the cholesky factorization of a symmetric positive-definite matrix
	 * \Output: _A := the lower triangular part L, so that $_A = L*L^T$
	 */
	MATRICE_HOST_INL static int spd(pointer _A, const size_type& _Size) {
		auto[M, N] = _Size;
#ifdef _DEBUG
		if (M != N) throw std::runtime_error("Non-sqaure matrix _A in _Lapack_kernel_impl<float>::spd(...).");
#endif // _DEBUG

#ifdef __use_mkl__
		return LAPACKE_spotrf(101, 'L', M, _A, N);
#else
		return flapk::_scholy(_A, N);
#endif
	}
	MATRICE_HOST_INL static int spd(const plview_type& _A) {
		return spd(std::get<2>(_A), { std::get<0>(_A), std::get<1>(_A) });
	}

	/**
	 * \computes the LU factorization of a general m-by-n matrix
	 * \Output: _A := L*U, _P := partial pivoting info
	 */
	MATRICE_HOST_INL static int lud(pointer _A, const size_type& _Size) {
		auto[M, N] = _Size;
		auto _P = new int[max(1, min(M, N))];

#ifdef __use_mkl__
		return LAPACKE_sgetrf(101, M, N, _A, max(1, N), _P);
#else
#ifdef _DEBUG
		if (M != N) throw std::runtime_error("Non-sqaure matrix _A in _Lapack_kernel_impl<float>::lud(...).");
#endif // _DEBUG
		return flak::_sLU(_A, _P, N)
#endif
	}
	MATRICE_HOST_INL static int lud(const plview_type& _A) {
		return lud(std::get<2>(_A), { std::get<0>(_A), std::get<1>(_A) });
	}

	/**
	 * \solves linear equations: A[N-by-N]*X[N-by-M] = B[N-by-M]
	 * \Output: $_B := \text(solutions) X$
	 */
	MATRICE_HOST_INL static int slv(const plview_type& _A, const plview_type& _B) {
#ifdef _DEBUG
		if (std::get<0>(_A) != std::get<1>(_A)) throw std::runtime_error("Non-sqaure matrix _A in _Lapack_kernel_impl<float>::slv(...).");
		if(std::get<1>(_A) != std::get<0>(_B)) throw std::runtime_error("Columns of _A .NEQ. rows of _B in _Lapack_kernel_impl<float>::slv(...).");
#endif
#ifdef __use_mkl__
		const auto N = std::get<0>(_A), M = std::get<1>(_B);
		auto _Ipiv = new int[max(1, N)];
		return LAPACKE_sgesv(101, N, M, std::get<2>(_A), N, _Ipiv, std::get<2>(_B), M);
#else
		throw std::runtime_error("Oops, no implementation is found.");
#endif
	}
};

/**
 *\Specialization for double-type.
 */
template<> struct _Lapack_kernel_impl<double> : _Lapack_kernel_impl_base<double> {
	/**
	 * \computes singular value decomposition
	 * \Output: $_A := U, _S := \Sigma, _Vt := V^T$
	 */
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

	/**
	 * \computes the cholesky factorization of a symmetric positive-definite matrix
	 * \Output: _A := the lower triangular part L, so that $_A = LL^T$
	 */
	MATRICE_HOST_INL static int spd(pointer _A, const size_type& _Size) {
		auto[M, N] = _Size;
#ifdef _DEBUG
		if (M != N) throw std::runtime_error("Non-sqaure matrix _A in _Lapack_kernel_impl<float>::spd(...).");
#endif // _DEBUG

#ifdef __use_mkl__
		return LAPACKE_dpotrf(101, 'L', M, _A, N);
#else
		return flapk::_dcholy(_A, N);
#endif
	}
	MATRICE_HOST_INL static int spd(const plview_type& _A) {
		return spd(std::get<2>(_A), { std::get<0>(_A), std::get<1>(_A) });
	}

	/**
	 * \computes the LU factorization of a general m-by-n matrix
	 * \Output: _A := L*U, _P := partial pivoting info
	 */
	MATRICE_HOST_INL static int lud(pointer _A, const size_type& _Size) {
		auto[M, N] = _Size;
		auto _P = new int[max(1, min(M, N))];

#ifdef __use_mkl__
		return LAPACKE_dgetrf(101, M, N, _A, max(1, N), _P);
#else
#ifdef _DEBUG
		if (M != N) throw std::runtime_error("Non-sqaure matrix _A in _Lapack_kernel_impl<float>::lud(...).");
#endif // _DEBUG
		return flak::_dLU(_A, _P, N)
#endif
	}
	MATRICE_HOST_INL static int lud(const plview_type& _A) {
		return lud(std::get<2>(_A), { std::get<0>(_A), std::get<1>(_A) });
	}

	/**
	 * \solves linear equations: A[N-by-N]*X[N-by-M] = B[N-by-M]
	 * \Output: $_B := \text(solutions) X$
	 */
	MATRICE_HOST_INL static int slv(const plview_type& _A, const plview_type& _B) {
#ifdef _DEBUG
		if (std::get<0>(_A) != std::get<1>(_A)) throw std::runtime_error("Non-sqaure matrix _A in _Lapack_kernel_impl<float>::slv(...).");
		if (std::get<1>(_A) != std::get<0>(_B)) throw std::runtime_error("Columns of _A .NEQ. rows of _B in _Lapack_kernel_impl<float>::slv(...).");
#endif
#ifdef __use_mkl__
		const auto N = std::get<0>(_A), M = std::get<1>(_B);
		auto _Ipiv = new int[max(1, N)];
		return LAPACKE_dgesv(101, N, M, std::get<2>(_A), N, _Ipiv, std::get<2>(_B), M);
#else
		throw std::runtime_error("Oops, no implementation is found.");
#endif
	}
};
_DETAIL_END DGE_MATRICE_END