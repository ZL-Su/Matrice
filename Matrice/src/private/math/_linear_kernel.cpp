/**********************************************************************
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
***********************************************************************/
#include <type_traits>
#include <core.hpp>
#include "private/math/_linear_kernel.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
namespace internal {
using _Size_t = long long;

/// <summary>
/// Solve linear system 'L*x=b', where 'L' a lower triangular matrix.
/// </summary>
/// <typeparam name="_Ptr">Pointer type</typeparam>
/// <param name="n">Number of rows or cols for the coeff. matrix.</param>
/// <param name="l">Pointer to n-by-n lower triangular matrix.</param>
/// <param name="y">Pointer to rhs vector, will be overwritten by the solution.</param>
/// <param name="stride">Number of columns of the rhs vector.</param>
template<typename _Ptr>
MATRICE_GLOBAL_INL void __tri_fwdsv_impl(_Size_t n, const _Ptr l, _Ptr y, _Size_t stride) noexcept {
	for (decltype(n) row_i = 0; row_i < n; ++row_i) { // over rows 0 to n
		auto& y_i = y[row_i * stride];
		const auto l_i = l + row_i * n;
		for (auto col_j = 0; col_j < row_i; ++col_j) { // over cols 0 to i
			y_i -= l_i[col_j] * y[col_j * stride];
		}
		y_i = safe_div(y_i, l_i[row_i]);
	}
}

template<typename _Ptr>
MATRICE_GLOBAL_INL void __tri_fwdsv_impl(_Size_t n, const _Ptr l, _Ptr y, int* p, _Size_t stride) noexcept {
	using value_type = primitive_type_t<_Ptr>;
	auto x = new value_type[n];
	for (auto i = 0; i < n; ++i) {
		x[i] = y[p[i] * stride];
	}
	for (auto i = 0; i < n; ++i) {
		y[i * stride] = x[i];
	}
	delete[] x;

	__tri_fwdsv_impl(n, l, y, stride);
}

/// <summary>
/// Solve linear system 'U*x=b', where 'U' an upper triangular matrix with unit diagonals.
/// </summary>
/// <typeparam name="_Ptr">Pointer type</typeparam>
/// <param name="n">Number of rows or cols for the coeff. matrix.</param>
/// <param name="u">Pointer to n-by-n upper triangular matrix.</param>
/// <param name="u">Pointer to rhs vector, will be overwritten by the solution.</param>
/// <param name="stride">Number of columns of the rhs vector.</param>
template<typename _Ptr>
MATRICE_GLOBAL_INL void __tri_bwdsv_impl(_Size_t n, const _Ptr u, _Ptr x, _Size_t stride) noexcept {
	for (decltype(n) i = n - 1; i >= 0; --i) { // over rows from n-1 to 0
		auto& x_i = x[i * stride];
		const auto u_i = u + i * n;
		for (auto j = i + 1; j < n; ++j) { // over cols from i+1 to n
			x_i -= u_i[j] * x[j * stride];
		}
	}
}

// \Perform Cholesky decomposition
template<typename _Ptr>
MATRICE_GLOBAL_INL int __spd_kernel_impl(_Ptr data, _Size_t n)noexcept {
	using value_type = primitive_type_t<_Ptr>;
	matrix_index_adapter idx{ n };

	int status = 1;
	auto _Sum = value_type(0);
	for (_Size_t r = 0; r < n; ++r) {
		for (_Size_t c = r; c < n; ++c) {
			_Sum = data[idx(r,c)];
			for (auto k = r - 1; k >= 0; --k) {
				_Sum -= data[idx(r, k)] * data[idx(c, k)];
			}
			if (r == c) {
				if (_Sum > value_type(0)) {
					data[idx(r, r)] = sqrt(_Sum);
				}
				else {
					return status = -r;
				}
			}
			else {
				data[idx(c, r)] = _Sum / data[idx(r, r)];
			}
		}
	}
	for (_Size_t r = 0; r < n; ++r)
		for (_Size_t c = 0; c < r; ++c)
			data[idx(c, r)] = value_type(0);

	return status;
}

/// <summary>
/// Compute inversion of a Cholesky factorized matrix.
/// </summary>
/// <typeparam name="_Ptr"></typeparam>
/// <param name="data">Input matrix</param>
/// <param name="inv">Output inversion</param>
/// <param name="n">Dimensinality.</param>
template<typename _Ptr>
MATRICE_GLOBAL_INL void __ispd_kernel_impl(const _Ptr data, _Ptr inv, _Size_t n) noexcept {
	using value_type = primitive_type_t<_Ptr>;
	value_type sum;
	for (auto r = 0; r < n; ++r) {
		for (auto c = 0; c <= r; ++c) {
			sum = r == c ? decltype(sum)(1) : 0;
			for (auto k = r - 1; k >= c; --k) {
				sum -= data[r * n + k] * inv[c * n + k];
			}
			inv[c * n + r] = sum / data[r * n + r];
		}
	}
	for (auto r = n - 1; r >= 0; --r) {
		for (auto c = 0; c <= r; ++c) {
			sum = (r < c ? 0. : inv[c * n + r]);
			for (auto k = r + 1; k < n; ++k)
				sum -= data[k * n + r] * inv[c * n + k];
			inv[r * n + c] = inv[c * n + r] = sum / data[r * n + r];
		}
	}
}

/// <summary>
/// Backward substitution to solve 'A*x = b' withe 'A' is symmetric position-definite.
/// </summary>
/// <typeparam name="_Ptr"></typeparam>
/// <param name="n">Dimensinality.</param>
/// <param name="lptr">Pointer to lower triangular matrix of 'A'.</param>
/// <param name="x">Pointer to instant column of the rhs vectors.</param>
/// <param name="stride">Number of columns of the rhs vectors.</param>
template<typename _Ptr>
MATRICE_GLOBAL_INL void __spd_bwdsv_impl(_Size_t n, const _Ptr lptr, _Ptr x, _Size_t stride) noexcept {
	// \solve: $L \cdot y = b$
	__tri_fwdsv_impl(n, lptr, x, stride);
	// \solve: $L^T \cdot x = y$
	for (auto r = n - 1; r >= 0; --r) {
		const auto x_idx = r * stride;
		auto t_sum = x[x_idx];
		for (auto c = r + 1; c < n; ++c) {
			t_sum -= lptr[c * n + r] * x[c * stride];
		}
		x[x_idx] = safe_div(t_sum, lptr[r * n + r]);
	}
}

// \Perform LU decomposition with pivoting
template<typename _Ptr>
MATRICE_GLOBAL int __lud_kernel_impl(_Size_t n, _Ptr data, int* piv) noexcept
{
	using value_type = primitive_type_t<_Ptr>;
	constexpr auto _Epsi{ std::numeric_limits<value_type>::epsilon() };
	matrix_index_adapter idx{ n };

	for (auto i = 0; i < n; ++i) {
		piv[i] = i;
	}
	auto& pivsign = piv[n] = { 1 };

	auto status{ 0 };
	auto pivot{ zero<value_type> };
	//loop over rows and columns...
	for (auto j = 0; j < n; j++) {
		auto p = j;
		pivot = data[idx(p, p)];
		for (auto i = j+1; i < n; i++) {
			if (auto tmp = abs(data[idx(i, j)]); tmp > pivot) {
				pivot = tmp; 
				p = i;
			}
		}
		if (abs(pivot) < sq(_Epsi)) {
			return status = -1;
		}
		if (j != p) {
			// interchange the pivot row and the current row
			for (auto col = 0; col < n; ++col) {
				swap(data[idx(j, col)], data[idx(p,col)]);
			}
			pivsign = -pivsign;
			swap(piv[j], piv[p]);
		}

		if (abs(data[idx(j, j)]) < _Epsi) data[idx(j, j)] = _Epsi;

		for (auto row = j + 1; row < n; ++row) {
			const auto tmp = data[idx(row,j)] /= data[idx(j,j)];
			for (auto col = j + 1; col < n; ++col)
				data[idx(row, col)] -= tmp * data[idx(j,col)];
		}
	}
	
	return status;
}

template<typename _Ptr>
MATRICE_GLOBAL int __svd_kernel_impl(_Size_t m, _Size_t n, _Ptr A, _Ptr W, _Ptr Vt) noexcept {
	using value_type = primitive_type_t<_Ptr>;
	int status = 0;

	return status;
}
}

template<>
MATRICE_GLOBAL int _Linear_spd_kernel(float* data, size_t n) noexcept {
	return internal::__spd_kernel_impl(data, n);
}
template<>
MATRICE_GLOBAL int _Linear_spd_kernel(double* data, size_t n) noexcept {
	return internal::__spd_kernel_impl(data, n);
}

template<>
MATRICE_GLOBAL void _Linear_ispd_kernel(float* data, float* inv, size_t n) noexcept {
	internal::__ispd_kernel_impl(data, inv, n);
}
template<>
MATRICE_GLOBAL void _Linear_ispd_kernel(double* data, double* inv, size_t n) noexcept {
	internal::__ispd_kernel_impl(data, inv, n);
}

template<> MATRICE_GLOBAL 
void _Linear_spd_bwd(size_t n, float* l, float* x, int stride) noexcept {
	internal::__spd_bwdsv_impl(n, l, x, stride);
}
template<> MATRICE_GLOBAL
void _Linear_spd_bwd(size_t n, double* l, double* x, int stride) noexcept {
	internal::__spd_bwdsv_impl(n, l, x, stride);
}

template<>
MATRICE_GLOBAL int _Linear_lud_kernel(size_t n, float* data, int* indx) noexcept {
	return internal::__lud_kernel_impl(n, data, indx);
}
template<>
MATRICE_GLOBAL int _Linear_lud_kernel(size_t n, double* data, int* indx) noexcept {
	return internal::__lud_kernel_impl(n, data, indx);
}

template<> MATRICE_GLOBAL
void _Linear_lud_sv(size_t n, float* lu, float* x, int* p, int stride) noexcept {
	internal::__tri_fwdsv_impl(n, lu, x, p, stride);
	internal::__tri_bwdsv_impl(n, lu, x, stride);
}
template<> MATRICE_GLOBAL
void _Linear_lud_sv(size_t n, double* lu, double* x, int*p, int stride) noexcept {
	internal::__tri_fwdsv_impl(n, lu, x, p, stride);
	internal::__tri_bwdsv_impl(n, lu, x, stride);
}
_DETAIL_END
DGE_MATRICE_END