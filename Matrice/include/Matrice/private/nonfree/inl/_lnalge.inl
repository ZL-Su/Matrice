/*  *************************************************************************
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
*	*************************************************************************/
#pragma once
#include <stdexcept>
#include "../_lnalge.h"
#include "../../math/_config.h"

DGE_MATRICE_BEGIN _DETAIL_BEGIN

template<ttag _Tag> struct transp_tag {};
template<> struct transp_tag<ttag::N> { 
#if   MATRICE_MATH_KERNEL == MATRICE_USE_MKL
	static constexpr auto value = CBLAS_TRANSPOSE::CblasNoTrans;
#else 
	static constexpr auto value = 111;
#endif
};
template<> struct transp_tag<ttag::Y> {
#if   MATRICE_MATH_KERNEL == MATRICE_USE_MKL
	static constexpr auto value = CBLAS_TRANSPOSE::CblasTrans;
#else 
	static constexpr auto value = 112;
#endif
};
template<ttag _Tag> MATRICE_HOST_INL constexpr auto transp_tag_v = transp_tag<_Tag>::value;

template<typename _Ty> struct _Blas_kernel_impl_base {
	using pointer = std::add_pointer_t<_Ty>;
	using size_type = tuple<int, int>;
	using plview_type = tuple<int, int, pointer, bool>;
#if   MATRICE_MATH_KERNEL == MATRICE_USE_MKL
	static constexpr auto layout = CBLAS_LAYOUT::CblasRowMajor;
#else 
	static constexpr auto layout = 101;
#endif
};
/**
 *\Specialization for float-type.
 */
template<> 
struct _Blas_kernel_impl<float> : _Blas_kernel_impl_base<float> {
	MATRICE_HOST_INL static auto dot(const pointer _x, const pointer _y, int _N, int _Incx = 1, int _Incy = 1) {
#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
		return cblas_sdot(_N, _x, _Incx, _y, _Incy);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports a kernel with preprocessor definition of MATRICE_MATH_KERNEL=MATRICE_USE_MKL.");
#endif
	}
	template<ttag Ta = ttag::N, ttag Tb = ttag::N>
	MATRICE_HOST_INL static auto mul(const pointer _A, const pointer _B, pointer _C, int _M, int _N, int _K) {
#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
		int _T = _K, lda = max(1,_K), ldb = max(1, _N);
		if constexpr (Ta == ttag::Y) {
			std::swap(_M, _K);
			lda = max(1, _M);
		}
		if constexpr (Tb == ttag::Y) {
			std::swap(_T, _N);
			ldb = max(1, _K);
		}
		cblas_sgemm(layout, transp_tag_v<Ta>, transp_tag_v<Tb>, _M, _N, _K, 1.f, _A, lda, _B, ldb, 0.f, _C, _N);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports a kernel with preprocessor definition of MATRICE_MATH_KERNEL=MATRICE_USE_MKL.");
#endif
	}
	template<ttag Ta = ttag::N, ttag Tb = ttag::N>
	MATRICE_HOST_INL static auto mul(const plview_type& _A, const plview_type& _B, plview_type& _C) {
		mul<Ta, Tb>(get<2>(_A), get<2>(_B), get<2>(_C), get<0>(_A), get<1>(_B), get<1>(_A));
	}

	MATRICE_HOST_INL static auto& gemm(const plview_type& _A, const plview_type& _B, plview_type _C) {
#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
		const auto[ma, na, a, ta] = _A;
		const auto[mb, nb, b, tb] = _B;
		auto c = get<2>(_C);
		if (mb == 1 || nb == 1) { //gemv
			const auto lda = max(1, na);
			cblas_sgemv(layout, 
				ta ? transp_tag_v<ttag::Y> : transp_tag_v<ttag::N>, 
				ma, na, 1.f, a, lda, b, 1, 1.f, c, 1);
		}
		else { //gemm
			int m = ma, k = nb, n = nb, t = k;
			int lda = max(1, k), ldb = max(1, n);
			if (ta) { std::swap(m, k); lda = max(1, m); }
			if (tb) { std::swap(t, n); ldb = max(1, k); }
			cblas_sgemm(layout,
				(ta ? transp_tag_v<ttag::Y> : transp_tag_v<ttag::N>),
				(tb ? transp_tag_v<ttag::Y> : transp_tag_v<ttag::N>),
				m, n, k, 1.f, a, lda, b, ldb, 0.f, c, n);
		}
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports a kernel with preprocessor definition of MATRICE_MATH_KERNEL=MATRICE_USE_MKL.");
#endif
		return (_C);
	}
};
/**
 *\Specialization for double-type.
 */
template<> struct _Blas_kernel_impl<double> : _Blas_kernel_impl_base<double> {
	MATRICE_HOST_INL static auto dot(const pointer _x, const pointer _y, int _N, int _Incx = 1, int _Incy = 1) {
#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
		return cblas_ddot(_N, _x, _Incx, _y, _Incy);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports a kernel with preprocessor definition of MATRICE_MATH_KERNEL=MATRICE_USE_MKL.");
#endif
	}
	template<ttag Ta = ttag::N, ttag Tb = ttag::N>
	MATRICE_HOST_INL static auto mul(const pointer _A, const pointer _B, pointer _C, int _M, int _N, int _K) {
#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
		int _T = _K, lda = max(1, _K), ldb = max(1, _N);
		if constexpr (Ta == ttag::Y) {
			std::swap(_M, _K);
			lda = max(1, _M);
	}
		if constexpr (Tb == ttag::Y) {
			std::swap(_T, _N);
			ldb = max(1, _K);
		}
		cblas_dgemm(layout, transp_tag_v<Ta>, transp_tag_v<Tb>, _M, _N, _K, 1., _A, _K, _B, _N, 0., _C, _N);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports a kernel with preprocessor definition of MATRICE_MATH_KERNEL=MATRICE_USE_MKL.");
#endif
	}
	template<ttag Ta = ttag::N, ttag Tb = ttag::N>
	MATRICE_HOST_INL static auto mul(const plview_type& _A, const plview_type& _B, plview_type& _C) {
		mul<Ta, Tb>(get<2>(_A), get<2>(_B), get<2>(_C), get<0>(_A), get<1>(_B), get<1>(_A));
	}

	MATRICE_HOST_INL static auto& gemm(const plview_type& _A, const plview_type& _B, plview_type _C) {
#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
		const auto[ma, na, a, ta] = _A;
		const auto[mb, nb, b, tb] = _B;
		auto c = get<2>(_C);
		if (mb == 1 || nb == 1) { //gemv
			const auto lda = max(1, na);
			cblas_dgemv(layout,
				ta ? transp_tag_v<ttag::Y> : transp_tag_v<ttag::N>,
				ma, na, 1.f, a, lda, b, 1, 1.f, c, 1);
		}
		else { //gemm
			int m = ma, k = nb, n = nb, t = k;
			int lda = max(1, k), ldb = max(1, n);
			if (ta) { std::swap(m, k); lda = max(1, m); }
			if (tb) { std::swap(t, n); ldb = max(1, k); }
			cblas_dgemm(layout,
				(ta ? transp_tag_v<ttag::Y> : transp_tag_v<ttag::N>),
				(tb ? transp_tag_v<ttag::Y> : transp_tag_v<ttag::N>),
				m, n, k, 1.f, a, lda, b, ldb, 0.f, c, n);
	}
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports a kernel with preprocessor definition of MATRICE_MATH_KERNEL=MATRICE_USE_MKL.");
#endif
		return (_C);
}
};

template<typename _Ty> struct _Lapack_kernel_impl_base {
	using pointer = std::add_pointer_t<_Ty>;
	using size_type = tuple<int, int>;
	using plview_type = tuple<int, int, pointer>;
#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
	static constexpr int layout = LAPACK_ROW_MAJOR;
#else
	static constexpr int layout = 101;
#endif
};
/**
 *\Specialization for float-type.
 */
template<> 
struct _Lapack_kernel_impl<float> : _Lapack_kernel_impl_base<float> {
	/**
	 * \computes singular value decomposition
	 * \Output: $_A := U, _S := \Sigma, _Vt := V^T$
	 */
	MATRICE_HOST_INL static int svd(pointer _A, pointer _S, pointer _Vt, const size_type& _Size) {
		const auto[M, N] = _Size;
#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
		float* _Superb = new float[(M < N ? M : N) - 1];
		return LAPACKE_sgesvd(layout, 'O', 'A',
			M, N, _A, N, _S, nullptr, 1, _Vt, N, _Superb);
#elif MATRICE_MATH_KERNEL == MATRICE_USE_FKL
		return flapk::_sgesvd(_A, _S, _Vt, M, N);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports two types of kernels with preprocessor definition of (1) MATRICE_MATH_KERNEL=MATRICE_USE_MKL, or (2) MATRICE_MATH_KERNEL=MATRICE_USE_FKL");
#endif
	}

	/**
	 * \computes the cholesky factorization of a symmetric positive-definite matrix
	 * \Output: _A := the lower triangular part L, so that $_A = L*L^T$
	 */
	MATRICE_HOST_INL static int spd(pointer _A, const size_type& _Size) {
		const auto[M, N] = _Size;
#ifdef _DEBUG
		DGELOM_CHECK(M == N, "Non-sqaure matrix _A in _Lapack_kernel_impl<float>::spd(...).");
#endif // _DEBUG

#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
		return LAPACKE_spotrf(layout, 'L', M, _A, N);
#elif MATRICE_MATH_KERNEL == MATRICE_USE_FKL
		return flapk::_scholy(_A, N);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports two types of kernels with preprocessor definition of (1) MATRICE_MATH_KERNEL=MATRICE_USE_MKL, or (2) MATRICE_MATH_KERNEL=MATRICE_USE_FKL");
#endif
	}
	MATRICE_HOST_INL static int spd(const plview_type& _A) {
		return spd(get<2>(_A), {get<0>(_A), get<1>(_A) });
	}
	template<typename _Mty, MATRICE_ENABLE_IF(is_matrix_v<_Mty>)>
	MATRICE_HOST_INL static int spd(const _Mty& _A) {
		return spd(_A.data(), _A.shape());
	}

	/**
	 * \computes the LU factorization of a general m-by-n matrix
	 * \Output: _A := L*U, _P := partial pivoting info
	 */
	MATRICE_HOST_INL static int lud(pointer _A, const size_type& _Size) {
		const auto[M, N] = _Size;
		auto _P = new int[max(1, min(M, N))];

#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
		return LAPACKE_sgetrf(layout, M, N, _A, max(1, N), _P);
#elif MATRICE_MATH_KERNEL == MATRICE_USE_FKL
#ifdef _DEBUG
		DGELOM_CHECK(M == N, "Non-sqaure matrix _A in _Lapack_kernel_impl<float>::lud(...).");
#endif // _DEBUG
		return flak::_sLU(_A, _P, N)
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports two types of kernels with preprocessor definition of (1) MATRICE_MATH_KERNEL=MATRICE_USE_MKL, or (2) MATRICE_MATH_KERNEL=MATRICE_USE_FKL");
#endif
	}
	MATRICE_HOST_INL static int lud(const plview_type& _A) {
		return lud(get<2>(_A), { get<0>(_A), get<1>(_A) });
	}

	/**
	 * \solves linear equations: A[N-by-N]*X[N-by-M] = B[N-by-M]
	 * \Output: $_B := \text(solutions) X$
	 */
	MATRICE_HOST_INL static int slv(const plview_type& _A, const plview_type& _B) {
#ifdef _DEBUG
		DGELOM_CHECK(get<0>(_A) == get<1>(_A),"Non-sqaure matrix _A in _Lapack_kernel_impl<float>::slv(...).");
		DGELOM_CHECK(get<1>(_A) == get<0>(_B),"Columns of _A .NEQ. rows of _B in _Lapack_kernel_impl<float>::slv(...).");
#endif
#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
		const auto N = get<0>(_A), M = get<1>(_B);
		auto _Ipiv = new int[max(1, N)];
		return LAPACKE_sgesv(layout, N, M, get<2>(_A), N, _Ipiv, get<2>(_B), M);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports a kernel with preprocessor definition of MATRICE_MATH_KERNEL=MATRICE_USE_MKL.");
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
		const auto[M, N] = _Size;
#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
		double* _Superb = new double[M < N ? M : N - 1];
		return LAPACKE_dgesvd(layout, 'O', 'A',
			M, N, _A, N, _S, nullptr, 1, _Vt, N, _Superb);
#elif MATRICE_MATH_KERNEL==MATRICE_USE_FKL
		return flapk::_dgesvd(_A, _S, _Vt, M, N);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports two types of kernels with preprocessor definition of (1) MATRICE_MATH_KERNEL=MATRICE_USE_MKL, or (2) MATRICE_MATH_KERNEL=MATRICE_USE_FKL");
#endif
	}

	/**
	 * \computes the cholesky factorization of a symmetric positive-definite matrix
	 * \Output: _A := the lower triangular part L, so that $_A = LL^T$
	 */
	MATRICE_HOST_INL static int spd(pointer _A, const size_type& _Size) {
		const auto[M, N] = _Size;
#ifdef _DEBUG
		DGELOM_CHECK(M == N,"Non-sqaure matrix _A in _Lapack_kernel_impl<float>::spd(...).");
#endif // _DEBUG

#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
		return LAPACKE_dpotrf(101, 'L', M, _A, N);
#elif MATRICE_MATH_KERNEL==MATRICE_USE_FKL
		return flapk::_dcholy(_A, N);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports two types of kernels with preprocessor definition of (1) MATRICE_MATH_KERNEL=MATRICE_USE_MKL, or (2) MATRICE_MATH_KERNEL=MATRICE_USE_FKL");
#endif
	}
	MATRICE_HOST_INL static int spd(const plview_type& _A) {
		return spd(get<2>(_A), { get<0>(_A), get<1>(_A) });
	}

	/**
	 * \computes the LU factorization of a general m-by-n matrix
	 * \Output: _A := L*U, _P := partial pivoting info
	 */
	MATRICE_HOST_INL static int lud(pointer _A, const size_type& _Size) {
		const auto[M, N] = _Size;
		auto _P = new int[max(1, min(M, N))];

#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
		return LAPACKE_dgetrf(layout, M, N, _A, max(1, N), _P);
#elif MATRICE_MATH_KERNEL==MATRICE_USE_FKL
#ifdef _DEBUG
		DGELOM_CHECK(M == N, "Non-sqaure matrix _A in _Lapack_kernel_impl<float>::lud(...).");
#endif // _DEBUG
		return flak::_dLU(_A, _P, N)
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports two types of kernels with preprocessor definition of (1) MATRICE_MATH_KERNEL=MATRICE_USE_MKL, or (2) MATRICE_MATH_KERNEL=MATRICE_USE_FKL");
#endif
	}
	MATRICE_HOST_INL static int lud(const plview_type& _A) {
		return lud(get<2>(_A), { get<0>(_A), get<1>(_A) });
	}

	/**
	 * \solves linear equations: A[N-by-N]*X[N-by-M] = B[N-by-M]
	 * \Output: $_B := \text(solutions) X$
	 */
	MATRICE_HOST_INL static int slv(const plview_type& _A, const plview_type& _B) {
#ifdef _DEBUG
		DGELOM_CHECK(get<0>(_A) != get<1>(_A), "Non-sqaure matrix _A in _Lapack_kernel_impl<float>::slv(...).");
		DGELOM_CHECK(get<1>(_A) != get<0>(_B), "Columns of _A .NEQ. rows of _B in _Lapack_kernel_impl<float>::slv(...).");
#endif
#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
		const auto N = get<0>(_A), M = get<1>(_B);
		auto _Ipiv = new int[max(1, N)];
		return LAPACKE_dgesv(layout, N, M, get<2>(_A), N, _Ipiv, get<2>(_B), M);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports a kernel with preprocessor definition of MATRICE_MATH_KERNEL=MATRICE_USE_MKL.");
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
		DGELOM_CHECK(M == N, "The coeff. _A in _Lapack_backward_impl<solver_type::CHD> must be a square matrix.");
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