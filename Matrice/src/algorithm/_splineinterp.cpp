/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
**************************************************************************/
#include "algs/interpolation/_splineinterp.h"

MATRICE_ALGS_BEGIN

#pragma region <!-- MACROS -->
#define _SPLCOEFF_CWISE_RECURSION(_A, _Z) \
for (size_t _Row = 1; _Row < _Height; ++_Row) { \
	_Buff.rview(_Row) = _Data.rview(_Row) + _Buff.rview(_Row-1)*_Z; \
} \
\
_Mycoeff.rview(_Height-1) = _A*(_Buff.rview(_Height-1)+_Buff.rview(_Height-2)*_Z); \
for (index_t _Row = _Height - 2; _Row >= 0; --_Row) { \
	_Mycoeff.rview(_Row) = _Z*(_Mycoeff.rview(_Row+1) - _Buff.rview(_Row)); \
}

#define _SPLCOEFF_RWISE_RECURSION(_A, _Z) \
for (size_t _Col = 1; _Col < _Width; ++_Col) { \
	_Buff.cview(_Col) = _Mycoeff.cview(_Col) + _Buff.cview(_Col - 1)*_Z; \
} \
\
_Mycoeff.cview(_Width-1) = _A*(_Buff.cview(_Width-1)+_Buff.cview(_Width-2)*_Z); \
for (index_t _Col = _Width - 2; _Col >= 0; --_Col) { \
	_Mycoeff.cview(_Col) = _Z*(_Mycoeff.cview(_Col+1) - _Buff.cview(_Col)); \
}
#pragma endregion

namespace internal {
	template<typename _Ty> struct _Itpar_base {
		static constexpr auto _Myeps = std::numeric_limits<_Ty>::epsilon();
	};
	template<typename _Ty, typename _Tag> struct _It_hypar {};

	template<typename _Ty> struct _It_hypar<_Ty, bicerp_tag> : _Itpar_base<_Ty> {
		MATRICE_HOST_INL static auto value() {
			constexpr auto _Myeps = _Itpar_base<_Ty>::_Myeps;

			const auto _Z = -2. + sqrt(3.);
			const auto _A = _Z / (_Z*_Z - 1);
			const auto _K = static_cast<index_t>(log(_Myeps) / log(abs(_Z))) + 2;

			return std::make_tuple(_Z, _A, _K);
		}
	};
	template<typename _Ty> struct _It_hypar<_Ty, biqerp_tag> : _Itpar_base<_Ty> {
		MATRICE_HOST_INL static auto value() {
			constexpr auto _Myeps = _Itpar_base<_Ty>::_Myeps;

			const auto _Z1 = -0.4305753470999737919, _Z2 = -0.04309628820326465382;
			const auto _A1 = _Z1 / (_Z1*_Z1 - 1), _A2 = _Z2 / (_Z2*_Z2 - 1);
			const auto _K1 = static_cast<index_t>(log(_Myeps)/log(abs(_Z1)))+2;
			const auto _K2 = static_cast<index_t>(log(_Myeps)/log(abs(_Z2)))+2;

			return std::make_tuple(_Z1, _Z2, _A1, _A2, _K1, _K2);
		}
	};
	template<typename _Ty> struct _It_hypar<_Ty, biserp_tag> : _Itpar_base<_Ty> {
		MATRICE_HOST_INL static auto value() {
			constexpr auto _Myeps = _Itpar_base<_Ty>::_Myeps;

			const auto _Z1 = -0.5352804307964381655;
			const auto _Z2 = -0.12255461519232669052;
			const auto _Z3 = -0.009148694809608276929;
			const auto _A1 = _Z1 / (_Z1*_Z1 - 1);
			const auto _A2 = _Z2 / (_Z2*_Z2 - 1);
			const auto _A3 = _Z3 / (_Z3*_Z3 - 1);
			const auto _K1 = static_cast<index_t>(log(_Myeps)/log(abs(_Z1)))+2;
			const auto _K2 = static_cast<index_t>(log(_Myeps)/log(abs(_Z2)))+2;
			const auto _K3 = static_cast<index_t>(log(_Myeps)/log(abs(_Z3)))+2;

			return std::make_tuple(_Z1, _Z2, _Z3, _A1, _A2, _A3, _K1, _K2, _K3);
		}
	};
}

template<typename _Ty> 
void _Spline_interpolation<_Ty, bicerp_tag>::_Coeff_impl() {
	const auto& _Data = *_Mybase::_Mydata;
	auto& _Mycoeff = _Mybase::_Mycoeff;

	auto[_Height, _Width] = _Data.shape();
	_Mycoeff.create(_Height, _Width, zero<value_type>);

	//initialization
	const auto[_Z, _A, _K] = internal::_It_hypar<value_type, _Mybase::category>::value();

	matrix_type _Buff(_Height, _Width);

	//Recursion over each column

	//pre-calculate the initial value for the first recursion
	auto _R0 = _Buff.rview(0) = zero<value_type>;
#pragma omp parallel if(_K > 100)
	{
#pragma omp for
		for (index_t _Row = 0; _Row < _K; ++_Row) {
			_R0 = _R0 + pow(_Z, _Row)*_Data.rview(_Row);
		}
	}
	_SPLCOEFF_CWISE_RECURSION(_A, _Z);

	//Recursion over each row

	//pre-calculate the initial value for the first recursion
	auto _C0 = _Buff.cview(0) = zero<value_type>;
#pragma omp parallel if(_K > 100)
	{
#pragma omp for
		for (index_t _Col = 0; _Col < _K; ++_Col) {
			_C0 = _C0 + _Mycoeff.cview(_Col) * pow(_Z, _Col);
		}
	}
	_SPLCOEFF_RWISE_RECURSION(_A, _Z);
}
template<typename _Ty> 
void _Spline_interpolation<_Ty, biqerp_tag>::_Coeff_impl() {
	const auto& _Data = *_Mybase::_Mydata;
	auto& _Mycoeff = _Mybase::_Mycoeff;

	const auto[_Height, _Width] = _Data.shape();
	_Mycoeff.create(_Height, _Width, zero<value_type>);

	//initialization
	const auto[_Z1, _Z2, _A1, _A2, _K1, _K2] = internal::_It_hypar<value_type, _Mybase::category>::value();

	matrix_type _Buff(_Height, _Width);

	// \Recursion over each column...

	//pre-calculate the initial value for the first recursion
	auto _R0 = _Buff.rview(0) = zero<value_type>;
#pragma omp parallel if(_K1 > 100)
	{
#pragma omp for
		for (index_t _Row = 0; _Row < _K1; ++_Row) {
			_R0 = _R0 + pow(_Z1, _Row)*_Data.rview(_Row);
		}
	}
	_SPLCOEFF_CWISE_RECURSION(_A1, _Z1);

	_R0 = zero<value_type>;
#pragma omp parallel if(_K2 > 100)
	{
#pragma omp for
		for (index_t _Row = 0; _Row < _K2; ++_Row) {
			_R0 = _R0 + pow(_Z2, _Row)*_Mycoeff.rview(_Row);
		}
	}
	_SPLCOEFF_CWISE_RECURSION(_A2, _Z2);

	// \Recursion over each row...

	//pre-calculate the initial value for the first recursion
	auto _C0 = _Buff.cview(0) = zero<value_type>;
#pragma omp parallel if(_K1 > 100)
	{
#pragma omp for
		for (index_t _Col = 0; _Col < _K1; ++_Col) {
			_C0 = _C0 + _Mycoeff.cview(_Col) * pow(_Z1, _Col);
		}
	}
	_SPLCOEFF_RWISE_RECURSION(_A1, _Z1);

	_C0 = zero<value_type>;
#pragma omp parallel if(_K2 > 100)
	{
#pragma omp for
		for (index_t _Col = 0; _Col < _K2; ++_Col) {
			_C0 = _C0 + _Mycoeff.cview(_Col) * pow(_Z2, _Col);
		}
	}
	_SPLCOEFF_RWISE_RECURSION(_A2, _Z2);
}
template<typename _Ty> 
void _Spline_interpolation<_Ty, biserp_tag>::_Coeff_impl() {
	const auto& _Data = *_Mybase::_Mydata;
	auto& _Mycoeff = _Mybase::_Mycoeff;

	const auto[_Height, _Width] = _Data.shape();
	_Mycoeff.create(_Height, _Width, zero<value_type>);

	//initialization
	const auto[_Z1, _Z2, _Z3, _A1, _A2, _A3, _K1, _K2, _K3] = internal::_It_hypar<value_type, _Mybase::category>::value();

	matrix_type _Buff(_Height, _Width);

	auto _R0 = _Buff.rview(0) = zero<value_type>;
#pragma omp parallel if(_K1 > 100)
	{
#pragma omp for
		for (index_t _Row = 0; _Row < _K1; ++_Row) {
			_R0 = _R0 + pow(_Z1, _Row)*_Data.rview(_Row);
		}
	}
	_SPLCOEFF_CWISE_RECURSION(_A1, _Z1);

	_R0 = zero<value_type>;
#pragma omp parallel if(_K2 > 100)
	{
#pragma omp for
		for (index_t _Row = 0; _Row < _K2; ++_Row) {
			_R0 = _R0 + pow(_Z2, _Row)*_Data.rview(_Row);
		}
	}
	_SPLCOEFF_CWISE_RECURSION(_A2, _Z2);

	_R0 = zero<value_type>;
#pragma omp parallel if(_K3 > 100)
	{
#pragma omp for
		for (index_t _Row = 0; _Row < _K3; ++_Row) {
			_R0 = _R0 + pow(_Z3, _Row)*_Data.rview(_Row);
		}
	}
	_SPLCOEFF_CWISE_RECURSION(_A3, _Z3);

	auto _C0 = _Buff.cview(0) = zero<value_type>;
#pragma omp parallel if(_K1 > 100)
	{
#pragma omp for
		for (index_t _Col = 0; _Col < _K1; ++_Col) {
			_C0 = _C0 + _Mycoeff.cview(_Col) * pow(_Z1, _Col);
		}
	}
	_SPLCOEFF_RWISE_RECURSION(_A1, _Z1);

	_C0 = zero<value_type>;
#pragma omp parallel if(_K2 > 100)
	{
#pragma omp for
		for (index_t _Col = 0; _Col < _K2; ++_Col) {
			_C0 = _C0 + _Mycoeff.cview(_Col) * pow(_Z2, _Col);
		}
	}
	_SPLCOEFF_RWISE_RECURSION(_A2, _Z2);

	_C0 = zero<value_type>;
#pragma omp parallel if(_K3 > 100)
	{
#pragma omp for
		for (index_t _Col = 0; _Col < _K3; ++_Col) {
			_C0 = _C0 + _Mycoeff.cview(_Col) * pow(_Z3, _Col);
		}
	}
	_SPLCOEFF_RWISE_RECURSION(_A3, _Z3);
}

template class _Spline_interpolation<float,  bicerp_tag>;
template class _Spline_interpolation<double, bicerp_tag>;
template class _Spline_interpolation<float,  biqerp_tag>;
template class _Spline_interpolation<double, biqerp_tag>;
template class _Spline_interpolation<float,  biserp_tag>;
template class _Spline_interpolation<double, biserp_tag>;

#undef _SPLCOEFF_CWISE_RECURSION
#undef _SPLCOEFF_RWISE_RECURSION

MATRICE_ALGS_END