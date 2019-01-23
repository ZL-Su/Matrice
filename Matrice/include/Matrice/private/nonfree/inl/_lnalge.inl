#pragma once
#include <stdexcept>
#include "../_lnalge.h"
#ifdef __use_mkl__
#include <mkl.h>
#else
#include <fkl.h>
#endif

DGE_MATRICE_BEGIN _DETAIL_BEGIN

template<ttag _Tag> struct transp_tag {};
template<> struct transp_tag<ttag::N> { 
	static constexpr auto value = CBLAS_TRANSPOSE::CblasNoTrans;
};
template<> struct transp_tag<ttag::Y> {
	static constexpr auto value = CBLAS_TRANSPOSE::CblasTrans;
};
template<ttag _Tag> MATRICE_HOST_INL constexpr auto transp_tag_v = transp_tag<_Tag>::value;

template<typename _Ty> struct _Blas_kernel_impl_base {
	using pointer = std::add_pointer_t<_Ty>;
	using size_type = std::tuple<int, int>;
	using plview_type = std::tuple<int, int, pointer>;
	static constexpr auto layout = CBLAS_LAYOUT::CblasRowMajor;
};
/**
 *\Specialization for float-type.
 */
template<> struct _Blas_kernel_impl<float> : _Blas_kernel_impl_base<float> {
	MATRICE_HOST_INL static auto dot(const pointer _x, const pointer _y, int _N, int _Incx = 1, int _Incy = 1) {
#ifdef __use_mkl__
		return cblas_sdot(_N, _x, _Incx, _y, _Incy);
#else
		throw std::runtime_error("Undefined kernel in _Blas_kernel_impl<float>::dot(...).");
#endif
	}
	template<ttag Ta = ttag::N, ttag Tb = ttag::N>
	MATRICE_HOST_INL static auto mul(const pointer _A, const pointer _B, pointer _C, int _M, int _N, int _K) {
#ifdef __use_mkl__
		cblas_sgemm(layout, transp_tag_v<Ta>, transp_tag_v<Tb>, _M, _N, _K, 1.f, _A, _K, _B, _N, 0.f, _C, _N);
#else
		throw std::runtime_error("Undefined kernel in _Blas_kernel_impl<float>::mul(...).");
#endif
	}
	template<ttag Ta = ttag::N, ttag Tb = ttag::N>
	MATRICE_HOST_INL static auto mul(const plview_type& _A, const plview_type& _B, const plview_type& _C) {
		mul<Ta, Tb>(std::get<2>(_A), std::get<2>(_B), std::get<2>(_C), std::get<0>(_A), std::get<1>(_B), std::get<1>(_A));
	}
};
/**
 *\Specialization for double-type.
 */
template<> struct _Blas_kernel_impl<double> : _Blas_kernel_impl_base<double> {
	MATRICE_HOST_INL static auto dot(const pointer _x, const pointer _y, int _N, int _Incx = 1, int _Incy = 1) {
#ifdef __use_mkl__
		return cblas_ddot(_N, _x, _Incx, _y, _Incy);
#else
		throw std::runtime_error("Undefined kernel in _Blas_kernel_impl<float>::mul(...).");
#endif
	}
	template<ttag Ta = ttag::N, ttag Tb = ttag::N>
	MATRICE_HOST_INL static auto mul(const pointer _A, const pointer _B, pointer _C, int _M, int _N, int _K) {
#ifdef __use_mkl__
		cblas_dgemm(layout, transp_tag_v<Ta>, transp_tag_v<Tb>, _M, _N, _K, 1., _A, _K, _B, _N, 0., _C, _N);
#else
		throw std::runtime_error("Undefined kernel in _Blas_kernel_impl<float>::mul(...).");
#endif
	}
	template<ttag Ta = ttag::N, ttag Tb = ttag::N>
	MATRICE_HOST_INL static auto mul(const plview_type& _A, const plview_type& _B, const plview_type& _C) {
		mul<Ta, Tb>(std::get<2>(_A), std::get<2>(_B), std::get<2>(_C), std::get<0>(_A), std::get<1>(_B), std::get<1>(_A));
	}
};

template<typename _Ty> struct _Lapack_kernel_impl_base {
	using pointer = std::add_pointer_t<_Ty>;
	using size_type = std::tuple<int, int>;
	using plview_type = std::tuple<int, int, pointer>;
	static constexpr int layout = LAPACK_ROW_MAJOR;
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
		return LAPACKE_sgesvd(layout, 'O', 'A',
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
		return LAPACKE_spotrf(layout, 'L', M, _A, N);
#else
		return flapk::_scholy(_A, N);
#endif
	}
	MATRICE_HOST_INL static int spd(const plview_type& _A) {
		return spd(std::get<2>(_A), { std::get<0>(_A), std::get<1>(_A) });
	}
	template<typename _Mty, typename = std::enable_if_t<is_matrix_v<_Mty>>>
	MATRICE_HOST_INL static int spd(const _Mty& _A) {
		return spd(_A.data(), _A.shape());
	}

	/**
	 * \computes the LU factorization of a general m-by-n matrix
	 * \Output: _A := L*U, _P := partial pivoting info
	 */
	MATRICE_HOST_INL static int lud(pointer _A, const size_type& _Size) {
		auto[M, N] = _Size;
		auto _P = new int[max(1, min(M, N))];

#ifdef __use_mkl__
		return LAPACKE_sgetrf(layout, M, N, _A, max(1, N), _P);
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
		return LAPACKE_sgesv(layout, N, M, std::get<2>(_A), N, _Ipiv, std::get<2>(_B), M);
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
		return LAPACKE_dgesvd(layout, 'O', 'A',
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
		return LAPACKE_dgetrf(layout, M, N, _A, max(1, N), _P);
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
		return LAPACKE_dgesv(layout, N, M, std::get<2>(_A), N, _Ipiv, std::get<2>(_B), M);
#else
		throw std::runtime_error("Oops, no implementation is found.");
#endif
	}
};

/**
 *\Specialization for the Chalosky factorization.
 */
template<> struct _Lapack_backward_impl<solver_type::CHD> {
	template<typename _Lhs, typename _Rhs>
	MATRICE_GLOBAL_INL static auto& eval(const _Lhs& _A, _Rhs& _X) {
		static_assert(is_matrix_v<_Lhs>, "_A in _Lapack_backward_impl<solver_type::CHD> must be a matrix type.");

		using iterator = typename _Lhs::iterator;
		using value_type = typename _Lhs::value_type;
		const auto M = _A.rows(), N = _A.cols();
		const auto NRhs = _X.cols();
#ifdef _DEBUG
		if (M != N) throw std::runtime_error("The coeff. _A in _Lapack_backward_impl<solver_type::CHD> must be a square matrix.");
#endif

		for (auto k = 0; k < NRhs; ++k) {
			auto _B = _X.cview(k);
			// \solve: L*y = b
			for (auto j = 0; j < M; ++j) {
				const auto _A_row = _A[j];
				auto _Sum_j = _B(j);
				for (auto i = 0; i < j; ++i) {
					_Sum_j -= _A_row[i] * _B(i);
				}
				_B(j) = _Sum_j / _A_row[j];
			}
			// \solve: U*x = y, where U = L^T
			for (auto j = M - 1; j >= 0; --j) {
				const auto _A_row = _A[j];
				auto _Sum_j = _B(j);
				for (auto i = j + 1; i < N; ++i) {
					_Sum_j -= _A_row[i] * _B(i);
				}
				_B(j) = _Sum_j / _A_row[j];
			}
		}
		return (_X);
	}
};
_DETAIL_END DGE_MATRICE_END