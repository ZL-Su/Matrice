#include <type_traits>
#include "math/_linear_fact.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ptr>
MATRICE_GLOBAL_INL int __spd_kernel(_Ptr data, size_t n) noexcept {
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
	return __spd_kernel(data, n);
}
template<>
MATRICE_GLOBAL int _Linear_spd_kernel(double* data, size_t n) noexcept {
	return __spd_kernel(data, n);
}

_DETAIL_END
DGE_MATRICE_END