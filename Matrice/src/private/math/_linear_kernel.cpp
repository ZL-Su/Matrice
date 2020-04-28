#include <type_traits>
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

	// \solve $u \cdot x = y$
	template<typename _Ptr>
	MATRICE_GLOBAL_INL void _tri_bwdsv_impl(_Size_t n, const _Ptr u, _Ptr x, _Size_t stride) noexcept {

	}

	// \Perform Cholesky decomposition
	template<typename _Ptr>
	MATRICE_GLOBAL_INL int __spd_kernel_impl(_Ptr data, _Size_t n)noexcept {
		using value_type = primitive_type_t<_Ptr>;
		int status = 1;
		auto _Sum = value_type(0);
		for (_Size_t r = 0; r < n; ++r) {
			for (_Size_t c = r; c < n; ++c) {
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
		for (_Size_t r = 0; r < n; ++r)
			for (_Size_t c = 0; c < r; ++c)
				data[c * n + r] = value_type(0);

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
	MATRICE_GLOBAL_INL int __lud_kernel_impl(_Ptr data, _Size_t n)noexcept {
		using value_type = primitive_type_t<_Ptr>;
		constexpr auto eps = std::numeric_limits<value_type>::epsilon();

		int status = 1;
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
	internal::__ispd_kernel_impl<float*>(data, inv, n);
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
_DETAIL_END
DGE_MATRICE_END