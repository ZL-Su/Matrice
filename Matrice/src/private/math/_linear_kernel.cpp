#include <type_traits>
#include "math/_linear_kernel.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ptr>
MATRICE_GLOBAL_INL int __spd_kernel_impl(_Ptr data, size_t n) noexcept {
	using size_type = long long;
	using value_type = primitive_type_t<_Ptr>;
	int status = 1;
	auto _Sum = value_type(0);
	for (size_type r = 0; r < n; ++r) {
		for (size_type c = r; c < n; ++c) {
			_Sum = data[r * n + c];
			for (auto k = r - 1; k >= 0; --k) {
				_Sum -= data[r * n + k] * data[c * n + k];
			}
			if (r == c) {
				if (_Sum > value_type(0)) {
					data[r * n + r] = sqrt(_Sum);
				}
				else {
					return status = -r;
				}
			}
			else {
				data[c * n + r] = _Sum / data[r * n + r];
			}
		}
	}
	for (size_type r = 0; r < n; ++r)
		for (size_type c = 0; c < r; ++c)
			data[c * n + r] = value_type(0);

	return status;
}
template<>
MATRICE_GLOBAL int _Linear_spd_kernel(float* data, size_t n) noexcept {
	return __spd_kernel_impl(data, n);
}
template<>
MATRICE_GLOBAL int _Linear_spd_kernel(double* data, size_t n) noexcept {
	return __spd_kernel_impl(data, n);
}

template<typename _Ptr>
MATRICE_GLOBAL_INL void __ispd_kernel_impl(const _Ptr data, _Ptr inv, size_t _n) noexcept {
	using size_type = long long;
	using value_type = primitive_type_t<_Ptr>;
	const size_type n = _n;
	value_type sum;
	for (auto r = 0; r < n; ++r) {
		for (auto c = 0; c <= r; ++c) {
			sum = r == c ? decltype(sum)(1) : 0;
			for (auto k = r - 1; k >= c; --k) {
				sum -= data[r*n+k] * inv[c*n+k];
			}
			inv[c*n+r] = sum / data[r*n+r];
		}
	}
	for (auto r = n - 1; r >= 0; --r) {
		for (auto c = 0; c <= r; ++c) {
			sum = (r < c ? 0. : inv[c*n+r]);
			for (auto k = r + 1; k < n; ++k)
				sum -= data[k*n+r] * inv[c*n+k];
			inv[r*n+c] = inv[c*n+r] = sum / data[r*n+r];
		}
	}
}
template<>
MATRICE_GLOBAL void _Linear_ispd_kernel(float* data, float* inv, size_t n) noexcept {
	__ispd_kernel_impl<float*>(data, inv, n);
}
template<>
MATRICE_GLOBAL void _Linear_ispd_kernel(double* data, double* inv, size_t n) noexcept {
	__ispd_kernel_impl(data, inv, n);
}

template<typename _Ptr>
MATRICE_GLOBAL_INL void __spd_bwd_kernel_impl(size_t n, const _Ptr lptr, _Ptr x, size_t stride) noexcept {
	// \solve: L*y = b
	for (auto r = 0; r < n; ++r) {
		const auto a_row = lptr + r * n;
		auto& t_sum = *(x + r * stride);
		for (auto c = 0; c < r; ++c) {
			t_sum -= a_row[c] * x[c * stride];
		}
		t_sum = safe_div(t_sum, a_row[r]);
	}
	// \solve: U*x = y, where U = L^T
	for (diff_t r = n - 1; r >= 0; --r) {
		const auto a_row = lptr + r * n;
		auto& t_sum = *(x + r * stride);
		for (auto c = r + 1; c < n; ++c) {
			t_sum -= a_row[c] * x[c * stride];
		}
		t_sum = safe_div(t_sum, a_row[r]);
	}
}
template<> MATRICE_GLOBAL 
void _Linear_spd_bwd(size_t n, float* data, float* x, int stride) noexcept {
	__spd_bwd_kernel_impl(n, data, x, stride);
}
template<> MATRICE_GLOBAL
void _Linear_spd_bwd(size_t n, double* data, double* x, int stride) noexcept {
	__spd_bwd_kernel_impl(n, data, x, stride);
}
_DETAIL_END
DGE_MATRICE_END