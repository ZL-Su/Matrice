#pragma once

#include "../_optim.h"
#include "../../../arch/ixpacket.h"

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

// \update on warp parameters
template<std::size_t _Order> struct _Compositional_warp_update{};

template<> struct _Compositional_warp_update<compile_time_size<>::val_1> {
	template<typename _Vecty>
	static MATRICE_GLOBAL_FINL auto& inv(const _Vecty& x, _Vecty y) {
		auto _Inv_det = 1 / ((1 + y[1])*(1 + y[5]) - y[2] * y[4]);
		auto _Val_0 = y[0] * (1 + y[5]) - y[3] * y[2];
		auto _Val_1 = y[3] * (1 + y[1]) - y[0] * y[4];
		auto _Val_2 = y[5] + 1, _Val_3 = y[4];
		auto _Val_4 = y[3] + 1, _Val_5 = y[2];

		y[0] = x[0] - ((1 + x[1])*_Val_0 - x[2] * _Val_1)*_Inv_det;
		y[1] = ((1 + x[1])*_Val_2 - x[2] * _Val_3)*_Inv_det - 1;
		y[2] = (x[2] * _Val_4 + (1 + x[1])*_Val_5)*_Inv_det;
		y[3] = x[3] - ((1 + x[5])*_Val_1 - x[4] * _Val_0)*_Inv_det;
		y[4] = (x[4] * _Val_2 - (1 + x[5])*_Val_3)*_Inv_det;
		y[5] = ((1 + x[5])*_Val_4 - x[4] * _Val_5)*_Inv_det - 1;

		return (y);
	}
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
	_Myref = _Diff * _Issd;

	if (_J.valid()) _J.get();

	_Myhess = _Myhess * (2 * _Issd);

	_Mysolver.forward();
}

template<typename _Derived> MATRICE_HOST_INL 
auto _Corr_optim_base<_Derived>::_Warp(const param_type& _Pars) {
#define _WITHIN_RANGE_OF_REFIMG(_OP) \
	if (x - _Bl >= 0 && x - _Bl >= 0 && y + _Bu < _Rows && x + _Bu < _Cols) \
		_OP; \
   else _Mycur(r,c) = zero_v<value_type>;

	using category = category_type_t<interp_type>;
	constexpr auto _Bl = _Corr_border_size<category>::lower;
	constexpr auto _Bu = _Corr_border_size<category>::upper;
	const auto [_Rows, _Cols] = _Mycur_itp().shape();
	const auto _Radius = static_cast<index_t>(_Myopt._Radius);
	auto _Cur_x = _Mypos.x + _Pars[0], _Cur_y = _Mypos.y + _Pars[3];

	auto _Mean = zero_v<value_type>;
	for (index_t j = -_Radius, r = 0; j <= _Radius; ++j, ++r) {
		auto dy = static_cast<value_type>(j);
		auto tx = _Cur_x + _Pars[2] * dy;
		auto ty = _Cur_y + _Pars[5] * dy + j;
		for (index_t i = -_Radius, c = 0; i <= _Radius; ++i, ++c) {
			auto dx = static_cast<value_type>(i);
			auto x = _Pars[1] * dx + tx + i, y = _Pars[4] * dx + ty;
			_WITHIN_RANGE_OF_REFIMG(_Mean += _Mycur(r, c) = _Mycur_itp(x,y));
		}
	}
	_Mean /= _Mycur.size();

	_Mycur = _Mycur - _Mean;
	auto _Issd = 1./sqrt((_Mycur*_Mycur).sum());
	_Mycur = _Mycur * _Issd;

	return (2 * _Issd);
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
			if constexpr (is_float32_v<value_type>) {
				using packet_t = simd::Packet_<value_type, 4>;
				auto a = packet_t{ dfdx, dfdx, dfdy, dfdy };
				auto b = (a * packet_t{ dx, dy, dx, dy }).begin();
				_Mybase::_Myjaco(j, i) = { dfdx, b[0], b[1], dfdy, b[2], b[3] };
			}
			else {
				_Mybase::_Myjaco(j, i) = {
				dfdx, dfdx*dx, dfdx*dy, dfdy, dfdy*dx, dfdy*dy
				};
			}
		}
	}

	return (_Mybase::_Myjaco);
}

template<typename _Ty, typename _Itag> MATRICE_HOST_INL
auto _Corr_invcomp_optim<_Ty, _Itag, 1>::_Update(param_type& _P) {
	// \warp
	const auto _Scal = _Mybase::_Warp(_P); 
	// \error image exp.
	auto _Diff_n = (_Mybase::_Myref - _Mybase::_Mycur)*_Scal; 
	// \steepest descent param. update
	param_type _Sdp = (_Diff_n*_Mybase::_Myjaco).reduce().t();
	// \solve warp param update
	_Sdp = -1.*_Mybase::_Mysolver.backward(_Sdp);

	auto _Error = (_Sdp*_Sdp).sum();
	if(_Error < 1.0E-3)
		// \inverse composition to update param.
		_P = _Compositional_warp_update<order>::inv(_P, _Sdp);

	return std::tuple((_Diff_n*_Diff_n).sum(), _Error);
}

_DETAIL_END } MATRICE_ALGS_END