/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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
#pragma once
#include <numeric>
#if (defined __enable_cuda__ && !defined __disable_cuda__)
#define __disable_simd__
#include <host_defines.h>
#define __global_inl__ __inline__ __host__ __device__
#else
#define __global_inl__ __forceinline
#endif

namespace dgelom { namespace privt {
template<typename _Ty> __global_inl__ int _inv2x2m(_Ty* A, _Ty* Ainv)
{
	double det = A[0] * A[3] - A[1] * A[2];

	if ((det < 0 ? -det : det) < 1.E-6) return 0;

	det = 1. / det; 
	volatile _Ty A0 = A[0];
	Ainv[0] =  A[3] * det, Ainv[1] = -A[1] * det;
	Ainv[2] = -A[2] * det, Ainv[3] =    A0 * det;

	return 1;
}

template<typename _Ty> __global_inl__ int _inv3x3m(_Ty* A, _Ty* Ainv)
{
	double det = A[0] * A[4] * A[8] + A[1] * A[5] * A[6] + A[2] * A[3] * A[7]
		- A[6] * A[4] * A[2] - A[7] * A[5] * A[0] - A[8] * A[3] * A[1];

	if ((det < 0 ? -det : det) < std::numeric_limits<double>::epsilon()) return 0;

	det = 1. / det;
	volatile _Ty A0 = A[0], A1 = A[1], A2 = A[2], A3 = A[3], A4 = A[4], A6 = A[6];
	Ainv[0] =  (A[4] * A[8] - A[7] * A[5]) * det;
	Ainv[1] = -(A[1] * A[8] - A[7] * A[2]) * det;
	Ainv[2] =  (A1   * A[5] - A[4] * A[2]) * det;

	Ainv[3] = -(A[3] * A[8] - A[6] * A[5]) * det;
	Ainv[4] =  (A0   * A[8] - A[6] * A2  ) * det;
	Ainv[5] = -(A0   * A[5] - A3   * A2  ) * det;

	Ainv[6] =  (A3 * A[7] - A[6] * A4) * det;
	Ainv[7] = -(A0 * A[7] - A6   * A1) * det;
	Ainv[8] =  (A0 * A4 -   A3   * A1) * det;

	return 1;
}
} }
