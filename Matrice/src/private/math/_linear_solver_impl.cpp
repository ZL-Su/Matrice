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
	auto [s, val] = detail::_Lapack_kernel_wrapper::gesvd(U, true);
	S = s;
	return solver_status{ val };
}
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<svd>(Matd_ref U, Matd_ref S, Matd_ref Vt) {
	auto [s, val] = detail::_Lapack_kernel_wrapper::gesvd(U, true);
	S = s;
	return solver_status{ val };
}
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<svd>(View_f32 U, View_f32 S, View_f32 Vt) {
	return solver_status();
}
template<> MATRICE_GLOBAL
solver_status _Lak_adapter<svd>(View_f64 U, View_f64 S, View_f64 Vt) {
	return solver_status();
}
template<> MATRICE_GLOBAL
solver_status _Inv_adapter<svd>(View_f32 U, View_f32 S, View_f32 V) {
	
	return solver_status();
}
template<> MATRICE_GLOBAL
solver_status _Inv_adapter<svd>(View_f64 U, View_f64 S, View_f64 V) {
	
	return solver_status();
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