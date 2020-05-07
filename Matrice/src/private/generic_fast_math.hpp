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

namespace dgelom { 
namespace internal {
template<typename _It, typename = std::enable_if_t<std::is_pointer_v<_It>>> 
__global_inl__ int _inv2x2m(const _It A, _It Ainv) {
	using value_type = typename std::pointer_traits<_It>::element_type;

	value_type _Det = A[0] * A[3] - A[1] * A[2];

	if ((_Det < 0 ? -_Det : _Det) < 1.E-6) return -1;

	_Det = 1. / _Det;
	Ainv[0] =  A[3] * _Det, Ainv[1] = -A[1] * _Det;
	Ainv[2] = -A[2] * _Det, Ainv[3] =  A[0] * _Det;

	return 0;
}

template<typename _It, typename = std::enable_if_t<std::is_pointer_v<_It>>>
__global_inl__ int _inv3x3m(const _It A, _It Ainv) {
	using value_type = typename std::pointer_traits<_It>::element_type;

	value_type _Det = A[0] * A[4] * A[8] + A[1] * A[5] * A[6] + A[2] * A[3] * A[7]
		- A[6] * A[4] * A[2] - A[7] * A[5] * A[0] - A[8] * A[3] * A[1];

	if ((_Det < 0 ? -_Det : _Det) < 1.E-6) return -1;

	_Det = 1. / _Det;
	volatile value_type A0 = A[0], A1 = A[1], A2 = A[2];
	volatile value_type A3 = A[3], A4 = A[4], A6 = A[6];

	Ainv[0] =  (A[4] * A[8] - A[7] * A[5]) * _Det;
	Ainv[1] = -(A[1] * A[8] - A[7] * A[2]) * _Det;
	Ainv[2] =  (A1   * A[5] - A[4] * A[2]) * _Det;

	Ainv[3] = -(A[3] * A[8] - A[6] * A[5]) * _Det;
	Ainv[4] =  (A0   * A[8] - A[6] * A2  ) * _Det;
	Ainv[5] = -(A0   * A[5] - A3   * A2  ) * _Det;

	Ainv[6] =  (A3 * A[7] - A[6] * A4) * _Det;
	Ainv[7] = -(A0 * A[7] - A6   * A1) * _Det;
	Ainv[8] =  (A0 * A4 -   A3   * A1) * _Det;

	return 0;
}
} }
