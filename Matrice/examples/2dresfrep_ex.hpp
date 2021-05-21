/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
**********************************************************************/
#pragma once

#include <algs/graphics.hpp>
#include <math/complex.h>

DGE_MATRICE_BEGIN
namespace example {
/// <summary>
/// \brief Generate 2D M-by-N sinusoidal signal with K frequencies.
/// </summary>
/// <param name="'nfreq'">Number of frequencies.</param>
/// <param name="'M'">Rows of the signal matrix to be generated.</param>
/// <param name="'N'">Columns of the signal matrix to be generated.</param>
/// <returns>A M-by-N signal matrix</returns>
auto _Gen_signal_2d(size_t nfreq, size_t M, size_t N = size_t()) noexcept {
	using value_t = float;

	N = N == size_t() ? M : N;

	const auto f1 = make_linspace(-0.5f, 0.5f, nfreq);
	const auto f2 = make_linspace(-0.45f, 0.45f, nfreq);
	const auto a = make_linspace(0.999f, 1.001f, nfreq);

	f1(0) = -0.4,   f1(1) = -0.35, f1(2) = -0.3,   f1(3) = -0.22;
	f1(4) = -0.16,  f1(5) = -0.04, f1(6) = 0.08,   f1(7) = 0.2;
	f2(0) = -0.269, f2(1) = -0.30, f2(2) = -0.277, f2(3) = -0.29;
	f2(4) = -0.26,  f2(5) = -0.25, f2(6) = -0.235, f2(7) = -0.29;
	a(0) = 1.039,   a(1) = 1.135,  a(2) = 0.222,   a(3) = 0.617;
	a(4) = 0.08,    a(5) = 0.959,  a(6) = 1.2,     a(7) = 0.217;

	complex_t<value_t, ::dynamic> Z(M, N);
	for (auto m = 0; m < M; ++m) {
		auto [Re, Im] = Z[m];
		for (const auto n : range<size_t>(0, N)) {
			auto& Reval = Re[n] = 0;
			auto& Imval = Im[n] = 0;
			for (const auto k : range<size_t>(0, nfreq)) {
				const auto x1 = 2 * pi<value_t>*f1(k) * m;
				const auto x2 = 2 * pi<value_t>*f2(k) * n;
				const auto c = complex_t<value_t>(sin(x1), cos(x1)) * complex_t<value_t>(sin(x2), cos(x2));
				const auto [real, imag] = a(k) * c;
				Reval += real, Imval += imag;
			}
		}
	}
	return forward<decltype(Z)>(Z);
}
}
DGE_MATRICE_END