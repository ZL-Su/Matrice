#pragma once

#include "../_optim.h"

MATRICE_ALGS_BEGIN _DETAIL_BEGIN namespace corr {

template<typename _Pty> MATRICE_HOST_INL 
auto _Patch_range(const _Pty& _Center, std::size_t _Radius) {
	int _L = floor<int>(_Center.x) - _Radius;
	int _R = floor<int>(_Center.x) + _Radius + 1;
	int _U = floor<int>(_Center.y) - _Radius;
	int _D = floor<int>(_Center.y) + _Radius + 1;

	return std::make_tuple(_L, _R, _U, _D);
}

template<typename _Derived> 
MATRICE_HOST_INL auto _Corr_optim_base<_Derived>::_Init() {
	const auto _Ksize = _Myopt._Radius << 1 | 1;
	_Myref.create(_Ksize, _Ksize);
	_Mycur.create(_Ksize, _Ksize);

	auto[_L, _R, _U, _D] = _Patch_range(_Mypos, _Myopt._Radius);
	for (auto y = _U; y < _D; ++y) {
		for (auto x = _L; x < _R; ++x) {
			_Myref(y - _U, x - _L) = _Myref_itp(x, y);
		}
	}
}

_DETAIL_END } MATRICE_ALGS_END