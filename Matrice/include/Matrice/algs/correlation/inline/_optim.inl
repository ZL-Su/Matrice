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
template<size_t _Order> struct _Compositional_warp_update{};

template<> 
struct _Compositional_warp_update<compile_time_size<>::_1> {
	template<typename _Vecty>
	static MATRICE_GLOBAL_FINL auto& inv(_Vecty& x, const _Vecty& y) {
		constexpr auto _One = one<typename _Vecty::value_t>;

		const auto y1p1 = _One + y[1], y5p1 = _One + y[5];
		const auto _Inv_of_det = _One / (y1p1*y5p1 - y[2]*y[4]);
		const auto a = y[2] * y[3] - y[0] * y5p1;
		const auto b = y[0] * y[4] - y[3] * y1p1;

		const auto x1p1 = _One + x[1], x5p1 = _One + x[5];
		const auto u  = _Inv_of_det*(x1p1*a    + x[2]*b) + x[0];
		const auto ux = _Inv_of_det*(x1p1*y5p1 - x[2]*y[4]) - _One;
		const auto uy = _Inv_of_det*(x[2]*y1p1 - y[2]*x1p1);
		const auto v  = _Inv_of_det*(x[4]*a    + x5p1*b) + x[3];
		const auto vx = _Inv_of_det*(x[4]*y5p1 - y[4]*x5p1);
		const auto vy = _Inv_of_det*(x5p1*y1p1 + x[4]*y[2]) - _One;

		x[0] = u, x[1] = ux, x[2] = uy, x[3] = v, x[4] = vx, x[5] = vy;

		return (x);
	}
};

// \base class implementation

template<typename _Derived> MATRICE_HOST_INL 
auto& _Corr_optim_base<_Derived>::_Init() {
	auto _J = std::async(std::launch::async, [&] {
		_MyJaco = static_cast<_Derived*>(this)->_Diff();
#if MATRICE_MATH_KERNEL==MATRICE_USE_NAT
		_Myhess = _MyJaco.t().mul(_MyJaco);
#else
		_Myhess = _MyJaco.inplace_mul<transp::Y>(_MyJaco);
#endif
	});

	const auto _Ksize = _Myopt._Radius << 1 | 1;
	_Myref.create(_Ksize, _Ksize);
	_Mycur.create(_Ksize, _Ksize);
	_Myerr.create(sqr(_Ksize), 1);

	const auto[_L, _R, _U, _D] = _Myopt.range<false>(_Mypos);
	auto _Mean = zero<value_type>;
	if (_Mypos.x == floor(_Mypos.x) && _Mypos.y == floor(_Mypos.y)) {
		_Myref = _Myref_itp.data().block(_L, _R, _U, _D).eval();
		_Mean = _Myref.sum() / _Myref.size();
	}
	else {
		for (auto y = _U; y < _D; ++y) {
			for (auto x = _L; x < _R; ++x) {
				_Mean += _Myref(y - _U, x - _L) = _Myref_itp(x, y);
			}
		}
	}

	auto _Diff = _Myref - _Mean;
	auto _Issd = one<value_type> / sqrt(sqr(_Diff).sum());
	_Myref = _Diff * _Issd;
	_Myissd = two<value_type>*_Issd;

	if (_J.valid()) _J.get();

	_Myhess = _Myhess*(_Myissd*_Issd);

	_Mysolver.forward();

	return (_Myref);
}

template<typename _Derived> MATRICE_HOST_INL 
auto& _Corr_optim_base<_Derived>::_Warp(const param_type& _Pars) {
#define _IF_WITHIN_RANGE(_OP) \
	if (x-_Bl>=0 && x-_Bl>=0 && y+_Bu<_Rows && x+_Bu<_Cols) \
		_OP; \
   else _Mycur(r,c) = zero<value_type>;

	using category = category_type_t<interp_type>;
	constexpr auto _Bl = _Corr_border_size<category>::lower;
	constexpr auto _Bu = _Corr_border_size<category>::upper;
	const auto [_Rows, _Cols] = _Mycur_itp().shape();
	const auto _Radius = static_cast<index_t>(_Myopt._Radius);
	const auto &u = _Pars[0],    &v = _Pars[3];
	const auto &dudx = _Pars[1], &dudy = _Pars[2];
	const auto &dvdx = _Pars[4], &dvdy = _Pars[5];
	const auto _Cur_x = _Mypos.x + u, _Cur_y = _Mypos.y + v;

	auto _Mean = zero<value_type>;
	const auto dvdyp1 = one<value_type> + _Pars[5];
	const auto dudxp1 = one<value_type> + _Pars[1];
	for (index_t j = -_Radius, r = 0; j <= _Radius; ++j, ++r) {
		const auto dy = static_cast<value_type>(j);
		const auto tx = _Cur_x + _Pars[2] * dy;
		const auto ty = _Cur_y + dvdyp1 * dy;
		for (index_t i = -_Radius, c = 0; i <= _Radius; ++i, ++c) {
			const auto dx = static_cast<value_type>(i);
			const auto x = dudxp1 * dx + tx, y = _Pars[4] * dx + ty;
			_IF_WITHIN_RANGE(_Mean += _Mycur(r, c) = _Mycur_itp(x,y));
		}
	}
	_Mean /= _Mycur.size();

	_Mycur = _Mycur - _Mean;
	const auto _Issd = 1./sqrt(sqr(_Mycur).sum());
	_Mycur = _Mycur * _Issd;

	return (_Mycur);
#undef _IF_WITHIN_RANGE
}

// \derived class implementation

template<typename _Ty, typename _Itag> MATRICE_HOST_INL
auto& _Corr_invcomp_optim<_Ty, _Itag, 1>::_Diff() {
	const auto _Size = _Mybase::_Myopt._Radius<<1|1;
	_Mybase::_MyJaco.create(sqr(_Size), _Mybase::DOF);

	auto[_L, _R, _U, _D] = _Mybase::_Myopt.range<true>(_Mybase::_Mypos);
	const auto _Off = -static_cast<value_type>(_Mybase::_Myopt._Radius);
	for (index_t iy = _U, j = 0; iy < _D; ++iy, ++j) {
		const auto dy = _Off + j, y = _Mybase::_Mypos.y + dy;
		for (index_t ix = _L, i = 0; ix < _R; ++ix, ++i) {
			const auto dx = _Off + i, x = _Mybase::_Mypos.x + dx;
			const auto[dfdx, dfdy] = _Mybase::_Myref_itp.grad({ x, y });

			auto p = _Mybase::_MyJaco[j * _Size + i];
			p[0] = dfdx, p[1] = dfdx * dx, p[2] = dfdx * dy;
			p[3] = dfdy, p[4] = dfdy * dx, p[5] = dfdy * dy;
		}
	}

	return (_Mybase::_MyJaco);
}

template<typename _Ty, typename _Itag> MATRICE_HOST_INL
auto _Corr_invcomp_optim<_Ty, _Itag, 1>::_Update(param_type& _P) {
	// \warp current subset
	_Mybase::_Mycur = _Mybase::_Warp(_P);

	// \error map
	_Mybase::_Myerr = (_Mybase::_Mycur-_Mybase::_Myref);
#ifdef _DEBUG
	typename _Mybase::matrix_type
		_Diff(_Mybase::_Myref.shape(), _Mybase::_Myerr.data());
#endif // _DEBUG

	// \steepest descent param. update
	param_type _Sdp = _Mybase::_MyJaco.t().mul(_Mybase::_Myerr);
	_Sdp = _Mybase::_Myissd*_Sdp;

	// \solve warp param update
	_Sdp = _Mybase::_Mysolver.backward(_Sdp);
	
	// \inverse composition to update param.
	_P = _Compositional_warp_update<order>::inv(_P, _Sdp);

	// \report least square correlation coeff. and param. error.
	return tuple(sqr(_Mybase::_Myerr).sum(), _Sdp.dot(_Sdp));
}

_DETAIL_END } MATRICE_ALGS_END