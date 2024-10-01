/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library for
Geometric Optical Sensing and Visual Intelligence.
Copyright(C) 2018-2024, Zhilong(Dgelom) Su, all rights reserved.

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
*********************************************************************/
#include "core/matrix.h"
#include "algs/spline/_spline.hpp"

DGE_MATRICE_BEGIN

#define MATRICE_BSPLCOEF_CWISE_RECURSION(_A, _Z) \
for (size_t _Row = 1; _Row < _H; ++_Row) { \
	_Buff.rview(_Row) = _Data.rview(_Row) + _Buff.rview(_Row-1)*_Z; \
} \
\
_Mycoef.rview(_H-1) = _A*(_Buff.rview(_H-1)+_Buff.rview(_H-2)*_Z); \
for (index_t _Row = _H - 2; _Row >= 0; --_Row) { \
	_Mycoef.rview(_Row) = _Z*(_Mycoef.rview(_Row+1)-_Buff.rview(_Row)); \
}

#define MATRICE_BSPLCOEF_RWISE_RECURSION(_A, _Z) \
for (size_t _Col = 1; _Col < _W; ++_Col) { \
	_Buff.cview(_Col) = _Mycoef.cview(_Col)+_Buff.cview(_Col-1)*_Z; \
} \
\
_Mycoef.cview(_W-1) = _A*(_Buff.cview(_W-1)+_Buff.cview(_W-2)*_Z); \
for (index_t _Col = _W - 2; _Col >= 0; --_Col) { \
	_Mycoef.cview(_Col) = _Z*(_Mycoef.cview(_Col+1)-_Buff.cview(_Col)); \
}

_DETAIL_BEGIN

template<typename _Ty> 
struct _Bspline_par_based {
	static constexpr auto _Myeps = std::numeric_limits<_Ty>::epsilon();
};

template<typename _Ty> 
struct _Bspline_par_3 : _Bspline_par_based<_Ty> {
	using _Bspline_par_based<_Ty>::_Myeps;

	MATRICE_HOST_INL static auto value() {
		MATRICE_USE_STD(make_tuple);

		constexpr auto _Z = -_Ty(2) + 1.732050808;
		constexpr auto _A = _Z / (sq(_Z) - 1);
		const auto _K = static_cast<index_t>(log(_Myeps)/log(abs(_Z)))+2;

		return make_tuple(_Z, _A, _K);
	}
};
template<typename _Ty> 
struct _Bspline_par_5 : _Bspline_par_based<_Ty> {
	using _Bspline_par_based<_Ty>::_Myeps;

	MATRICE_HOST_INL static auto value() {
		MATRICE_USE_STD(make_tuple);

		constexpr auto _Z1 = -0.4305753470999737919;
		constexpr auto _Z2 = -0.04309628820326465382;
		constexpr auto _A1 = _Z1 / (_Z1 * _Z1 - 1);
		constexpr auto _A2 = _Z2 / (_Z2 * _Z2 - 1);
		const auto _K1 = static_cast<index_t>(log(_Myeps)/log(abs(_Z1)))+2;
		const auto _K2 = static_cast<index_t>(log(_Myeps)/log(abs(_Z2)))+2;

		return make_tuple(_Z1, _Z2, _A1, _A2, _K1, _K2);
	}
};
template<typename _Ty> 
struct _Bspline_par_7 : _Bspline_par_based<_Ty> {
	using _Bspline_par_based<_Ty>::_Myeps;

	MATRICE_HOST_INL static auto value() {
		MATRICE_USE_STD(make_tuple);

		constexpr auto _Z1 = -0.5352804307964381655;
		constexpr auto _Z2 = -0.12255461519232669052;
		constexpr auto _Z3 = -0.009148694809608276929;
		constexpr auto _A1 = _Z1 / (_Z1 * _Z1 - 1);
		constexpr auto _A2 = _Z2 / (_Z2 * _Z2 - 1);
		constexpr auto _A3 = _Z3 / (_Z3 * _Z3 - 1);
		const auto _K1 = static_cast<index_t>(log(_Myeps)/log(abs(_Z1)))+2;
		const auto _K2 = static_cast<index_t>(log(_Myeps)/log(abs(_Z2)))+2;
		const auto _K3 = static_cast<index_t>(log(_Myeps)/log(abs(_Z3)))+2;

		return make_tuple(_Z1, _Z2, _Z3, _A1, _A2, _A3, _K1, _K2, _K3);
	}
};

/// <summary>
/// \brief Bicubic B-spline constructor
/// </summary>
template<typename _Ty>
_Bspline<_Ty, bicubic_tag>::_Bspline(bool _Zp) {
	_Mybase::_Mycoef = make_shared_matrix<_Ty>();
	_Mybase::_Myzp = _Zp;
}
template<typename _Ty>
_Bspline<_Ty, bicubic_tag>::_Bspline(const matrix_type& _Data, bool _Zp) {
	_Mybase::_Myzp = _Zp;
	_Mybase::_Mycoef = make_shared_matrix<_Ty>(_Data.shape());
	_Precompute(_Data);
}

template<typename _Ty>
void _Bspline<_Ty, bicubic_tag>::set(const matrix_type& _Data) {
	// require zero padding
	if (_Data.shape() != _Mybase::_Mycoef->shape()) {
		constexpr auto _Bound= ~-(4 >> 1) + -~(4 >> 1);
		const auto [_H, _W, _] = _Data.shape();
		_Mybase::_Mycoef->create(_H, _W);
	}
	_Precompute(_Data);
}

/// <summary>
/// \brief Precompute coefficient matrix
/// </summary>
template<typename _Ty>
void _Bspline<_Ty, bicubic_tag>::_Precompute(const matrix_type& _Data) {
#ifdef MATRICE_DEBUG
	DGELOM_CHECK(_Mybase::_Mycoef, 
		"Uninitialized coefficient matrix in _Bspline<_Ty, bicubic_tag>.");
#endif
	decltype(auto) _Mycoef = *_Mybase::_Mycoef;
	const auto [_Z, _A, _K] = _Bspline_par_3<value_type>::value();
	const auto [_H, _W, _] = _Mycoef.shape();
	matrix_type _Buff(_H, _W);

	// Precompute initial values for the first recursion
	auto _R0 = _Buff.rview(0) = value_type(0);
#pragma omp parallel if(_K > 5)
	{
#pragma omp for
		for (index_t _Row = 0; _Row < _K; ++_Row) {
			_R0 = _R0 + pow(_Z, _Row) * _Data.rview(_Row);
		}
	}

	// Recursion over each column
	MATRICE_BSPLCOEF_CWISE_RECURSION(_A, _Z);

	// Precompute initial values for the first recursion
	auto _C0 = _Buff.cview(0) = value_type(0);
#pragma omp parallel if(_K > 5)
	{
#pragma omp for
		for (index_t _Col = 0; _Col < _K; ++_Col) {
			_C0 = _C0 + _Mycoef.cview(_Col) * pow(_Z, _Col);
		}
	}

	// Recursion over each row
	MATRICE_BSPLCOEF_RWISE_RECURSION(_A, _Z);
}

template class _Bspline<float, bicubic_tag>;
template class _Bspline<double, bicubic_tag>;
_DETAIL_END

#undef MATRICE_BSPLCOEF_CWISE_RECURSION
#undef MATRICE_BSPLCOEF_RWISE_RECURSION

DGE_MATRICE_END