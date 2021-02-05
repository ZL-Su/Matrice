/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
#include "math/_linear_solver.hpp"
#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
#include "math/kernel_wrapper.hpp"
#else
#include "math/_linear_kernel.hpp"
#endif
#include "math/_linear_kernel.hpp" // for debug

#include "core/matrix.h"

DGE_MATRICE_BEGIN
_INTERNAL_BEGIN
using Matf_ref = Matrix<float>&;
using Matd_ref = Matrix<double>&;
using View_f32 = detail::_Matrix_block<float>;
using View_f64 = detail::_Matrix_block<double>;
using View_int = detail::_Matrix_block<int>;

// specialization for svd_op

/// <summary>
/// Perform Singular Value Decomposition (SVD).
/// </summary>
/// <param name="'U'">Source data matrix, overwritten with the U-matrix.</param>
/// <param name="'S'">Vector holds the singular values in descending order.</param>
/// <param name="'Vt'">Transpose of the V-matrix.</param>
/// <returns> '0' for success </returns>
template<> MATRICE_GLOBAL 
solver_status _Lak_adapter<svd>(Matf_ref U, Matf_ref S, Matf_ref Vt) {
	const auto lyt = U.allocator().fmt();
	const auto ldu = min(U.rows(), U.cols());
	remove_all_t<Matf_ref> supb(ldu);
	int ret = internal::_lapack_gesvd(lyt, 'S', 'S', U.rows(), U.cols(), U.data(), U.cols(), 
		S.data(), U.data(), ldu, Vt.data(), Vt.cols(), supb.data());
	return solver_status{ ret };
}
/// <summary>
/// Perform Singular Value Decomposition (SVD).
/// </summary>
/// <param name="'U'">Source data matrix, overwritten with the U-matrix.</param>
/// <param name="'S'">Vector holds the singular values in descending order.</param>
/// <param name="'Vt'">Transpose of the V-matrix.</param>
/// <returns> '0' for success </returns>
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<svd>(Matd_ref U, Matd_ref S, Matd_ref Vt) {
	const auto lyt = U.allocator().fmt();
	const auto ldu = min(U.rows(), U.cols());
	remove_all_t<Matd_ref> supb(ldu);
	int ret = internal::_lapack_gesvd(lyt, 'S', 'S', U.rows(), U.cols(), U.data(), U.cols(), 
		S.data(), U.data(), ldu, Vt.data(), Vt.cols(), supb.data());
	return solver_status{ ret };
}
/// <summary>
/// Perform Singular Value Decomposition (SVD).
/// </summary>
/// <param name="'U'">Source data matrix, overwritten with the U-matrix.</param>
/// <param name="'S'">Vector holds the singular values in descending order.</param>
/// <param name="'Vt'">Transpose of the V-matrix.</param>
/// <returns> '0' for success </returns>
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<svd>(View_f32 U, View_f32 S, View_f32 Vt) {
	const auto ldu = min(U.rows(), U.cols());
	matrix_f32 supb(ldu);
	int ret = internal::_lapack_gesvd(101, 'S', 'S', U.rows(), U.cols(), U[0], U.cols(), 
		S[0], U[0], ldu, Vt[0], Vt.cols(), supb.data());
	return solver_status{ ret };
}
/// <summary>
/// Perform Singular Value Decomposition (SVD).
/// </summary>
/// <param name="'U'">Source data matrix, overwritten with the U-matrix.</param>
/// <param name="'S'">Vector holds the singular values in descending order.</param>
/// <param name="'Vt'">Transpose of the V-matrix.</param>
/// <returns> '0' for success </returns>
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<svd>(View_f64 U, View_f64 S, View_f64 Vt) {
	const auto ldu = min(U.rows(), U.cols());
	matrix_f64 supb(ldu);
	int ret = internal::_lapack_gesvd(101, 'S', 'S', U.rows(), U.cols(), U[0], U.cols(), 
		S[0], U[0], ldu, Vt[0], Vt.cols(), supb.data());
	return solver_status{ ret };
}
template<> MATRICE_GLOBAL solver_status 
_Bwd_adapter<svd>(View_f32 U, View_f32 S, View_f32 Vt, View_f32 b, View_f32 x){
	auto tmp = transpose(U).mul(b).eval();
	const auto _Thresh = decltype(tmp)::eps * sqrt<float>(U.size()) * S(0) / 2;
	for (auto n = 0; n < x.size(); ++n) {
		tmp(n) = safe_div(tmp(n), S(n), _Thresh);
	}
	x = transpose(Vt).mul(tmp);
	return solver_status{ 0 };
}
template<> MATRICE_GLOBAL solver_status 
_Bwd_adapter<svd>(View_f64 U, View_f64 S, View_f64 Vt, View_f64 b, View_f64 x) {
	auto tmp = transpose(U).mul(b).eval();
	const auto _Thresh = decltype(tmp)::eps * sqrt<double>(U.size()) * S(0) / 2;
	for (auto n = 0; n < x.size(); ++n) {
		tmp(n) = safe_div(tmp(n), S(n), _Thresh);
	}
	x = transpose(Vt).mul(tmp);
	return solver_status{ 0 };
}
template<> MATRICE_GLOBAL
solver_status _Inv_adapter<svd>(View_f32 U, View_f32 S, View_f32 Vt, View_f32 Inv) {
	const auto _Thresh = std::numeric_limits<float>::epsilon()*sqrt<float>(U.size())*S(0) / 2;
	for (auto r = 0; r < U.rows(); ++r) {
		const auto pU = U[r];
		for (auto c = 0; c < U.cols(); ++c) {
			auto sum = View_f32::value_t(0);
			for (auto k = 0; k < Vt.rows(); ++k) {
				sum += pU[k] * safe_div(Vt[k][c], S(k), _Thresh);
			}
			Inv[c][r] = sum;
		}
	}
	return solver_status{ 0 };
}
template<> MATRICE_GLOBAL
solver_status _Inv_adapter<svd>(View_f64 U, View_f64 S, View_f64 Vt, View_f64 Inv) {
	const auto _Thresh = std::numeric_limits<double>::epsilon()*sqrt<double>(U.size())*S(0) / 2;
	for (auto r = 0; r < U.rows(); ++r) {
		const auto pU = U[r];
		for (auto c = 0; c < U.cols(); ++c) {
			auto sum = View_f64::value_t(0);
			for (auto k = 0; k < Vt.rows(); ++k) {
				sum += pU[k] * safe_div(Vt[k][c], S(k), _Thresh);
			}
			Inv[c][r] = sum;
		}
	}
	return solver_status{ 0 };
}

// specialization for spt_op
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<spt>(Matf_ref A) {
	solver_status status;
	status.value = detail::_Linear_spd_kernel(A.data(), A.rows());
	return status;
}
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<spt>(Matd_ref A) {
	solver_status status;
	status.value = detail::_Linear_spd_kernel(A.data(), A.rows());
	return status;
}
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<spt>(View_f32 A) {
	solver_status status;
	status.value = detail::_Linear_spd_kernel(A[0], A.rows());
	return status;
}
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<spt>(View_f64 A) {
	solver_status status;
	status.value = detail::_Linear_spd_kernel(A[0], A.rows());
	return status;
}
template<> MATRICE_GLOBAL
solver_status _Bwd_adapter<spt>(View_f32 L, View_f32 B) {
	const auto _Stride = B.cols();
	for (auto _Off = 0; _Off < _Stride; ++_Off) {
		auto b = B[0] + _Off;
		detail::_Linear_spd_bwd(L.rows(), L[0], b, _Stride);
	}
	return solver_status{ 0 };
}
template<> MATRICE_GLOBAL
solver_status _Bwd_adapter<spt>(View_f64 L, View_f64 B) {
	const auto _Stride = B.cols();
	for (auto _Off = 0; _Off < _Stride; ++_Off) {
		auto b = B[0] + _Off;
		detail::_Linear_spd_bwd(L.rows(), L[0], b, _Stride);
	}
	return solver_status{ 0 };
}
template<> MATRICE_GLOBAL
solver_status _Inv_adapter<spt>(View_f32 L, View_f32 INV) {
	detail::_Linear_ispd_kernel(L[0], INV[0], L.rows());
	return solver_status{ 0 };
}
template<> MATRICE_GLOBAL
solver_status _Inv_adapter<spt>(View_f64 L, View_f64 INV) {
	detail::_Linear_ispd_kernel(L[0], INV[0], L.rows());
	return solver_status{ 0 };
}

template<> MATRICE_GLOBAL
solver_status _Imp_adapter(View_f32 A, View_f32 b, View_f32 x) {
	b = matmul(A, x) - b;
	return solver_status{0 };
}
template<> MATRICE_GLOBAL
solver_status _Imp_adapter(View_f64 A, View_f64 b, View_f64 x) {
	b = matmul(A, x) - b;
	return solver_status{ 0 };
}

/// <summary>
/// Computes the LU factorization of a general m-by-n matrix.
/// </summary>
/// <param name="'A'">View to matrix A.</param>
/// <param name="'piv'">View to permutation vector.</param>
/// <returns>Status with value '0' for sucess.</returns>
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<lud>(View_f32 A, View_int piv) {
	return solver_status{
#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
		LAPACKE_sgetrf(MKL_ROW_MAJOR, 
		A.rows(), A.cols(), A.data(), A.cols(), piv.data())
#else
		detail::_Linear_lud_kernel(A.rows(), A.data(), piv.data())
#endif
	};
}
/// <summary>
/// Computes the LU factorization of a general m-by-n matrix.
/// </summary>
/// <param name="'A'">View to matrix A.</param>
/// <param name="'piv'">View to permutation vector.</param>
/// <returns>Status with value '0' for sucess.</returns>
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<lud>(View_f64 A, View_int piv) {
	return solver_status{
#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
		LAPACKE_dgetrf(MKL_ROW_MAJOR, 
		A.rows(), A.cols(), A.data(), A.cols(), piv.data())
#else
		detail::_Linear_lud_kernel(A.rows(), A.data(), piv.data())
#endif
	};
}

/// <summary>
/// Solve LU factorized linear system 'LU*X=B' with data type "float".
/// </summary>
/// <param name="'LU'">LU-factorized coeff. matrix.</param>
/// <param name="'B'">Right-hand side vector, with multi-columns optionally.</param>
/// <param name="'p'">Permutation vector.</param>
/// <returns>'solver_status.value = 0'</returns>
template<> MATRICE_GLOBAL
solver_status _Bwd_adapter<lud>(View_f32 LU, View_f32 B, View_int p) {
#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
	return solver_status{
		LAPACKE_sgetrs(MKL_ROW_MAJOR, 'N', LU.rows(), B.cols(), LU.data(),
		LU.cols(), p.data(), B.data(), B.cols())
	};
#else
	const auto _Stride = B.cols();
	for (auto _Off = 0; _Off < _Stride; ++_Off) {
		auto b = B[0] + _Off;
		detail::_Linear_lud_sv(LU.rows(), LU[0], b, p[0], _Stride);
	}
	return solver_status{ 0 };
#endif
}

/// <summary>
/// Solve LU factorized linear system 'LU*X=B' with data type "double".
/// </summary>
/// <param name="'LU'">LU-factorized coeff. matrix.</param>
/// <param name="'B'">Right-hand side vector, with multi-columns optionally.</param>
/// <param name="'p'">Permutation vector.</param>
/// <returns>'solver_status.value = 0'</returns>
template<> MATRICE_GLOBAL
solver_status _Bwd_adapter<lud>(View_f64 LU, View_f64 B, View_int p) {
#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
	return solver_status{
		LAPACKE_dgetrs(MKL_ROW_MAJOR, 'N', LU.rows(), B.cols(), LU.data(),
		LU.cols(), p.data(), B.data(), B.cols())
	};
#else
	const auto _Stride = B.cols();
	for (auto _Off = 0; _Off < _Stride; ++_Off) {
		auto b = B[0] + _Off;
		detail::_Linear_lud_sv(LU.rows(), LU[0], b, p[0], _Stride);
	}
	return solver_status{ 0 };
#endif
}
_INTERNAL_END
DGE_MATRICE_END