/***********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2022, Zhilong(Dgelom) Su, all rights reserved.

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

#include <core/matrix.h>

MATRICE_ALG_BEGIN(optim)
/// <summary>
/// \brief CLASS TEMPLATE LAD
/// Least Absolute Deviations fitting via ADMM, which solves the following problem:
///                             minimize ||Ax - b||_1
/// </summary>
/// <typeparam name="_Ty">Scalar type</typeparam>
/// <typeparam name="_My">Matrix type</typeparam>
template<typename _Ty, typename _My>
class LAD {
public:
	using value_type = _Ty;
	using matrix_t = _My;
	struct options_type {
		// \brief Over-relaxation parameter with typical values between 1.0 and 1.8.
		value_type alpha = 1.4;

		// \brief Augmented Lagrangian parameter.
		value_type rho = 1.0;

		// \brief Max iteration number.
		size_t max_iter = 1000;

		value_type abs_tol = 1.E-4;
		value_type rel_tol = 1.E-2;
	};
	using options_t = options_type;

	explicit LAD(const matrix_t& A, const matrix_t& b, const options_t& opt = {})
		:_Myopt{ opt }, _Mya(A), _Myb(b) {
	}

	MATRICE_HOST_INL auto solve() const {
		const auto m = _Mya.rows();
		const auto n = _Mya.cols();

		auto x = matrix_t::zeros(n, 1);
		auto z = matrix_t::zeros(m, 1);
		auto u = matrix_t::zeros(m, 1);
		for (auto k = 0; k < _Myopt.max_iter; ++k) {
			if (k > 1) {

			}
			else {

			}
		}

		return x;
	}

private:
	matrix_t _Mya, _Myb;
	options_t _Myopt;
};
MATRICE_ALG_END(optim)