/***********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#pragma once
#include <util/_macros.h>

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
/// <summary>
/// A fast version for native matrix multiplication.
/// </summary>
/// <typeparam name="T"> scalar value type </typeparam>
/// <param name="C">: Matrix to be computed</param>
/// <param name="A">: Left-hand side matrix</param>
/// <param name="B">: Right-hand side matrix</param>
/// <param name="M">: Row dimension for C and A</param>
/// <param name="N">: Column dimensions for C and B</param>
/// <param name="K">: Column dimension for A, row dimension for B</param>
/// <returns>Pointer to the resulted matrix C</returns>
template<typename T>
MATRICE_GLOBAL_INL T* _Fast_native_matmul(T* C, const T* A, const T* B, size_t M, size_t N, size_t K) noexcept {
	for (auto i = 0; i < M; ++i) {
		const auto A_i = A + i * K;
		auto C_i = C + i * N;
		for (auto k = 0; k < K; ++k) {
			const auto A_ik = A_i[k];
			const auto B_k = B + k * N;
			for (auto j = 0; j < N; ++j) {
				C_i[j] = A_ik * B_k[j];
			}
		}
	}
	return C;
}
_DETAIL_END
DGE_MATRICE_END