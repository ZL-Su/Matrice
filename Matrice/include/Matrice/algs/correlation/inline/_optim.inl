#pragma once

#include "../_optim.h"

MATRICE_ALGS_BEGIN _DETAIL_BEGIN namespace corr {

// \retrieve border size for each interpolation alg.
template<typename _Itptag> struct _Corr_border_size {};
template<> struct _Corr_border_size<_TAG bicspl_tag> {
	static constexpr auto lower = 1, upper = 2;
};
template<> struct _Corr_border_size<_TAG biqspl_tag> {
	static constexpr auto lower = 2, upper = 3;
};
template<> struct _Corr_border_size<_TAG bisspl_tag> {
	static constexpr auto lower = 3, upper = 4;
};

// \base class implementation

template<typename _Derived> MATRICE_HOST_INL 
auto _Corr_optim_base<_Derived>::_Init() {
	auto _J = std::async(std::launch::async, [&] {
		static_cast<_Derived*>(this)->_Diff();
		_Myhess = _Myjaco.t().mul(_Myjaco).reduce();
	});

	const auto _Ksize = _Myopt._Radius << 1 | 1;
	_Myref.create(_Ksize, _Ksize);
	_Mycur.create(_Ksize, _Ksize);

	auto[_L, _R, _U, _D] = _Myopt.range<false>(_Mypos);
	for (auto y = _U; y < _D; ++y) {
		for (auto x = _L; x < _R; ++x) {
			_Myref(y - _U, x - _L) = _Myref_itp(x, y);
		}
	}

	auto _Mean = _Myref.sum() / _Myref.size();
	auto _Diff = _Myref - _Mean;
	auto _Issd = 1 / sqrt((_Diff*_Diff).sum());
	_Myref = (_Myref - _Mean)*_Issd;

	if (_J.valid()) _J.get();

	_Myhess = _Myhess * (2 * _Issd);

	_Mysolver.forward();
}

template<typename _Derived> MATRICE_HOST_INL 
auto _Corr_optim_base<_Derived>::_Warp(const param_type& _Pars) {

}

// \derived class implementation

template<typename _Ty, typename _Itag> MATRICE_HOST_INL
auto& _Corr_invcomp_optim<_Ty, _Itag, 1>::_Diff() {
	const auto _Size = _Mybase::_Myopt._Radius<<1|1;
	_Mybase::_Myjaco.create(_Size, _Size);

	auto[_L, _R, _U, _D] = _Mybase::_Myopt.range<true>(_Mybase::_Mypos);
	auto _Off = -static_cast<value_type>(_Mybase::_Myopt._Radius);
	for (index_t iy = _U, j = 0; iy < _D; ++iy, ++j) {
		auto dy = _Off + j, y = _Mybase::_Mypos.y + dy;
		for (index_t ix = _L, i = 0; ix < _R; ++ix, ++i) {
			auto dx = _Off + i, x = _Mybase::_Mypos.x + dx;

			auto[dfdx, dfdy] = _Mybase::_Myref_itp.grad({ x, y });
			_Mybase::_Myjaco(j, i) = {
				dfdx, dfdx*dx, dfdx*dy, dfdy, dfdy*dx, dfdy*dy
			};
		}
	}

	return (_Mybase::_Myjaco);
}

template<typename _Ty, typename _Itag> MATRICE_HOST_INL
auto& _Corr_invcomp_optim<_Ty, _Itag, 1>::_Update(param_type& _P) {



	return (_P);
}

_DETAIL_END } MATRICE_ALGS_END