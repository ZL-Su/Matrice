/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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
#include "math/_linear_kernel.hpp"

#include "core/matrix.h"

DGE_MATRICE_BEGIN
_INTERNAL_BEGIN
using Matf_ref = Matrix<float>&;
using Matd_ref = Matrix<double>&;
using View_f32 = detail::_Matrix_block<float>;
using View_f64 = detail::_Matrix_block<double>;

// specialization for svd_op
template<> MATRICE_GLOBAL 
solver_status _Lak_adapter<svd>(Matf_ref U, Matf_ref S, Matf_ref Vt) {
	const auto lyt = U.allocator().fmt();
	const auto ldu = min(U.rows(), U.cols());
	remove_all_t<Matf_ref> supb(ldu);
	auto ret = internal::_lapack_gesvd(lyt, 'S', 'S', U.rows(), U.cols(), U.data(), U.cols(), S.data(), U.data(), ldu, Vt.data(), Vt.cols(), supb.data());
	return solver_status{ ret };
}
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<svd>(Matd_ref U, Matd_ref S, Matd_ref Vt) {
	const auto lyt = U.allocator().fmt();
	const auto ldu = min(U.rows(), U.cols());
	remove_all_t<Matd_ref> supb(ldu);
	auto ret = internal::_lapack_gesvd(lyt, 'S', 'S', U.rows(), U.cols(), U.data(), U.cols(), S.data(), U.data(), ldu, Vt.data(), Vt.cols(), supb.data());
	return solver_status{ ret };
}
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<svd>(View_f32 U, View_f32 S, View_f32 Vt) {
	const auto ldu = min(U.rows(), U.cols());
	matrix_f32 supb(ldu);
	auto ret = internal::_lapack_gesvd(101, 'S', 'S', U.rows(), U.cols(), U[0], U.cols(), S[0], U[0], ldu, Vt[0], Vt.cols(), supb.data());
	return solver_status{ ret };
}
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<svd>(View_f64 U, View_f64 S, View_f64 Vt) {
	const auto ldu = min(U.rows(), U.cols());
	matrix_f64 supb(ldu);
	auto ret = internal::_lapack_gesvd(101, 'S', 'S', U.rows(), U.cols(), U[0], U.cols(), S[0], U[0], ldu, Vt[0], Vt.cols(), supb.data());
	return solver_status{ ret };
}
template<> MATRICE_GLOBAL
solver_status _Bwd_adapter<svd>(View_f32 U, View_f32 S, View_f32 Vt, View_f32 b, View_f32 x){
	auto tmp = transpose(U).mul(b).eval();
	const auto _Thresh = decltype(tmp)::eps * sqrt<float>(U.size()) * S(0) / 2;
	for (auto n = 0; n < x.size(); ++n) {
		tmp(n) = safe_div(tmp(n), S(n), _Thresh);
	}
	x = transpose(Vt).mul(tmp);
	return solver_status{ 1 };
}
template<> MATRICE_GLOBAL
solver_status _Bwd_adapter<svd>(View_f64 U, View_f64 S, View_f64 Vt, View_f64 b, View_f64 x) {
	auto tmp = transpose(U).mul(b).eval();
	const auto _Thresh = decltype(tmp)::eps * sqrt<double>(U.size()) * S(0) / 2;
	for (auto n = 0; n < x.size(); ++n) {
		tmp(n) = safe_div(tmp(n), S(n), _Thresh);
	}
	x = transpose(Vt).mul(tmp);
	return solver_status{ 1 };
}
template<> MATRICE_GLOBAL
solver_status _Inv_adapter<svd>(View_f32 U, View_f32 S, View_f32 Vt, View_f32 Inv) {
	auto V = transpose(Vt).eval();
	const auto _Thresh = decltype(V)::eps * sqrt<float>(U.size()) * S(0) / 2;
	V.rview(0) = V.rview(0) * safe_div(1, S(0), _Thresh);
	V.rview(1) = V.rview(1) * safe_div(1, S(1), _Thresh);
	V.rview(2) = V.rview(2) * safe_div(1, S(2), _Thresh);
	Inv = V.mul(transpose(U));
	return solver_status{ 1 };
}
template<> MATRICE_GLOBAL
solver_status _Inv_adapter<svd>(View_f64 U, View_f64 S, View_f64 Vt, View_f64 Inv) {
	auto V = transpose(Vt).eval();
	const auto _Thresh = decltype(V)::eps * sqrt<double>(U.size()) * S(0) / 2;
	V.rview(0) = V.rview(0) * safe_div(1, S(0), _Thresh);
	V.rview(1) = V.rview(1) * safe_div(1, S(1), _Thresh);
	V.rview(2) = V.rview(2) * safe_div(1, S(2), _Thresh);
	Inv = V.mul(transpose(U));
	return solver_status{ 1 };
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
	return solver_status();
}
template<> MATRICE_GLOBAL
solver_status _Bwd_adapter<spt>(View_f64 L, View_f64 B) {
	const auto _Stride = B.cols();
	for (auto _Off = 0; _Off < _Stride; ++_Off) {
		auto b = B[0] + _Off;
		detail::_Linear_spd_bwd(L.rows(), L[0], b, _Stride);
	}
	return solver_status();
}
template<> MATRICE_GLOBAL
solver_status _Inv_adapter<spt>(View_f32 L, View_f32 INV) {
	detail::_Linear_ispd_kernel(L[0], INV[0], L.rows());
	return solver_status();
}
template<> MATRICE_GLOBAL
solver_status _Inv_adapter<spt>(View_f64 L, View_f64 INV) {
	detail::_Linear_ispd_kernel(L[0], INV[0], L.rows());
	return solver_status();
}
_INTERNAL_END
DGE_MATRICE_END