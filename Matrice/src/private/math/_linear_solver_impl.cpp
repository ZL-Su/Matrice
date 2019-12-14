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
#include "math/_linear_fact.hpp"
#endif

#include "core/matrix.h"
#include <iostream>

DGE_MATRICE_BEGIN
_INTERNAL_BEGIN
using Matf_ref = Matrix<float>&;
using Matd_ref = Matrix<double>&;

// specialization for svd_op
template<> MATRICE_GLOBAL 
void _Lak_adapter<svd>(Matf_ref U, Matf_ref S, Matf_ref Vt) {
	std::cout << "specialization for svd_op\n";
}
template<> MATRICE_GLOBAL
void _Lak_adapter<svd>(Matd_ref U, Matd_ref S, Matd_ref Vt) {
	
}

// specialization for spt_op
template<> MATRICE_GLOBAL
void _Lak_adapter<spt>(Matf_ref A) {
	auto sts = detail::_Linear_spd_kernel(A.data(), A.rows());
}
template<> MATRICE_GLOBAL
void _Lak_adapter<spt>(Matd_ref A) {
	auto sts = detail::_Linear_spd_kernel(A.data(), A.rows());
}
_INTERNAL_END
DGE_MATRICE_END