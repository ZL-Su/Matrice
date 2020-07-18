#include <type_traits>
#include <core.hpp>
#include "math/_linear_kernel.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
namespace internal {
using _Size_t = long long;
// \solve $L \cdot y = b$
template<typename _Ptr>
MATRICE_GLOBAL_INL void _tri_fwdsv_impl(_Size_t n, const _Ptr l, _Ptr y, _Size_t stride) noexcept {
	y[0] = safe_div(y[0], l[0][0]);
	for (_Size_t i = 1; i < n; ++i) {
		auto& y_i = y[i * stride];
		const auto l_i = l + i * n;
		for (auto j = 0; j < i; ++j) {
			y_i -= l_i[j] * y[j * stride];
		}
		y_i = safe_div(y_i, l_i[i]);
	}
}

// \solve $U \cdot x = y$
template<typename _Ptr>
MATRICE_GLOBAL_INL void _tri_bwdsv_impl(_Size_t n, const _Ptr u, _Ptr x, _Size_t stride) noexcept {

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

template<typename _Ptr>
MATRICE_GLOBAL_INL void __spd_bwdsv_impl(_Size_t n, const _Ptr lptr, _Ptr x, _Size_t stride) noexcept {
	// \solve: $L \cdot y = b$
	for (auto r = 0; r < n; ++r) {
		const auto a_row = lptr + r * n;
		const auto x_idx = r * stride;
		auto t_sum = x[x_idx];
		for (auto c = 0; c < r; ++c) {
			t_sum -= a_row[c] * x[c * stride];
		}
		x[x_idx] = safe_div(t_sum, a_row[r]);
	}
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
MATRICE_GLOBAL int __lud_kernel_impl(_Size_t n, _Ptr data, int* indx)noexcept {
	using value_type = primitive_type_t<_Ptr>;
	constexpr auto _Epsi = std::numeric_limits<value_type>::epsilon();
	matrix_index_adapter idx{ n };

	for (auto row = 0; row < n; ++row) indx[row] = row;
	indx[n] = 1;

	int status = 0;
	for (auto k = 0; k < n; ++k) {
		auto piv = zero<value_type>;
		decltype(n) major_row = k;
		for (auto row = k; row < n; ++row) {
			auto tmp = abs(data[idx(row,k)]);
			if (tmp > piv) {
				piv = tmp; 
				major_row = row;
			}
			if (abs(piv) < sq(_Epsi)) status = -row;
		}
		if (k != major_row) {
			// interchange the major row and the k-th row
			for (auto col = 0; col < n; ++col) {
				auto tmp = data[major_row*n + col];
				data[major_row * n + col] = data[idx(k, col)];
				data[idx(k, col)] = tmp;
			}
			std::swap(indx[major_row], indx[k]);
			indx[n] = -indx[n];
		}

		if (abs(data[idx(k, k)]) < _Epsi) data[idx(k, k)] = _Epsi;

		for (auto row = k + 1; row < n; ++row) {
			auto tmp = data[idx(row,k)] /= data[idx(k,k)];
			for (auto col = k + 1; col < n; ++col)
				data[idx(row, col)] -= tmp * data[idx(k,col)];
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
void _Linear_spd_bwd(size_t n, float* data, float* x, int stride) noexcept {
	internal::__spd_bwdsv_impl(n, data, x, stride);
}
template<> MATRICE_GLOBAL
void _Linear_spd_bwd(size_t n, double* data, double* x, int stride) noexcept {
	internal::__spd_bwdsv_impl(n, data, x, stride);
}

template<>
MATRICE_GLOBAL int _Linear_lud_kernel(size_t n, float* data, int* indx) noexcept {
	return internal::__lud_kernel_impl(n, data, indx);
}
template<>
MATRICE_GLOBAL int _Linear_lud_kernel(size_t n, double* data, int* indx) noexcept {
	return internal::__lud_kernel_impl(n, data, indx);
}
_DETAIL_END
DGE_MATRICE_END