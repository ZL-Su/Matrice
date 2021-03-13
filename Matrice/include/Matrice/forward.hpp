/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for
more detail.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once
#include "util/_std_wrapper.h"

namespace dgelom {
namespace detail {
	template<typename _Ty, int _M, int _N> class Matrix_;
	template<typename _Ty, int _Dim> class Vec_;
	template<typename _Ty> class Vec3_;
	template<typename _Ty> class Vec4_;
}

template<typename _Ty, int _RowsAtCompileTime, int _ColsAtCompileTime=_RowsAtCompileTime> 
using Matrix_ = detail::Matrix_<_Ty, _RowsAtCompileTime, _ColsAtCompileTime>;

//\matrix type with host dynamic memory allocator
template<typename _Ty> 
using Matrix = Matrix_<_Ty, 0>;

template<typename _Ty> class Scalar;

#ifdef MATRICE_ENABLE_CUDA
//\matrix type with unified memory allocator
template<typename T, size_t _Options = rmaj | gene> 
using Umatrix = Matrix_<T,
	compile_time_size<>::RunTimeDeducedOnDevice,
	compile_time_size<>::RunTimeDeducedOnHost>;

//\matrix type with device memory allocator
template<typename T, size_t _Options = rmaj | gene> 
using Dmatrix = Matrix_<T,
	compile_time_size<>::RunTimeDeducedOnDevice,
	compile_time_size<>::RunTimeDeducedOnDevice>;
#endif

//\dynamic matrix type with single(32)/double(64) floating point type
using matrix_f32 = Matrix<float>;
using matrix_f64 = Matrix<double>;

//\N-dimensional array with host managed memory
template<typename T, int _N> 
using array_n = Matrix_<T, _N, 1>;

template<typename _Ty, int _Size>
using Array_ = Matrix_<_Ty, _Size, (_Size > 0 ? 1 : _Size)> ;

//\brief shared matrix pointer
template<typename _Ty>
using shared_matrix_t = shared_ptr<Matrix<_Ty>>;
}