#include "../../include/Matrice/algs/interpolation/_splineinterp.h"

MATRICE_ALGS_BEGIN
#pragma region <!-- These codes will be deprecated -->
template<typename _Ty, size_t _Options>
void BicubicSplineInterp<_Ty, _Options>::_Bspline_coeff(const matrix_t& _Data) {
	const auto _Width = _Data.cols(), _Height = _Data.rows();
	m_coeff.create(_Height, _Width, zero<value_t>::value);

	//initialization
	const value_t _Z = -2. + sqrt(value_t(3));
	const value_t _A = _Z / (_Z*_Z - 1);
	const index_t _K0 = static_cast<index_t>(log(base_t::m_eps) / log(abs(_Z))) + 1;

	matrix_t _Buff(_Height, _Width);

	//Recursion over each column

	//pre-calculate the initial value for the first recursion
	matrix_t _R0(1, _Width, zero<value_t>::value);
#pragma omp parallel if(_K0 > 100)
{
#pragma omp for
	for (index_t _Row = 0; _Row < _K0; ++_Row) {
		_R0 = _R0 + pow(_Z, _Row)*_Data.rview(_Row);
	}
}

	_Buff.rview(0) = _R0;
	for (size_t _Row = 1; _Row < _Height; ++_Row) {
		_Buff.rview(_Row) = _Data.rview(_Row) + _Buff.rview(_Row - 1)*_Z;
	}

	m_coeff.rview(_Height - 1) = _A * (_Buff.rview(_Height - 1) + _Buff.rview(_Height - 2)*_Z);
	for (index_t _Row = _Height - 2; _Row >= 0; --_Row) {
		m_coeff.rview(_Row) = _Z * (m_coeff.rview(_Row + 1) - _Buff.rview(_Row));
	}

	//Recursion over each row

	//pre-calculate the initial value for the first recursion
	matrix_t _C0(_Height, 1, zero<value_t>::value);
#pragma omp parallel if(_K0 > 100)
{
#pragma omp for
	for (index_t _Col = 0; _Col < _K0; ++_Col) {
		_C0 = _C0 + m_coeff.cview(_Col) * pow(_Z, _Col);
	}
}

	_Buff.cview(0) = _C0;
	for (size_t _Col = 1; _Col < _Width; ++_Col) {
		_Buff.cview(_Col) = m_coeff.cview(_Col) + _Buff.cview(_Col - 1)*_Z;
	}

	m_coeff.cview(_Width - 1) = _A * (_Buff.cview(_Width - 1) + _Buff.cview(_Width - 2)*_Z);
	for (index_t _Col = _Width - 2; _Col >= 0; --_Col) {
		m_coeff.cview(_Col) = _Z*(m_coeff.cview(_Col + 1) - _Buff.cview(_Col));
	}
}
template class BicubicSplineInterp<float, INTERP | BSPLINE | BICUBIC>;
template class BicubicSplineInterp<double, INTERP | BSPLINE | BICUBIC>;
#pragma endregion

template<typename _Ty>
void _Spline_interpolation<_Ty, _BICBSPL>::_Coeff_impl() {
	const auto& _Data = _Mybase::_Mydata;
	auto& _Mycoeff = _Mybase::_Mycoeff;

	auto[_Width, _Height] = _Data.shape();
	_Mycoeff.create(_Height, _Width, zero<value_type>::value);

	//initialization
	const value_type _Z = -2. + sqrt(value_type(3));
	const value_type _A = _Z / (_Z*_Z - 1);
	const index_t _K0 = static_cast<index_t>(log(_Mybase::_Myeps) / log(abs(_Z))) + 1;

	matrix_type _Buff(_Height, _Width);

	//Recursion over each column

	//pre-calculate the initial value for the first recursion
	matrix_type _R0(1, _Width, zero<value_type>::value);
#pragma omp parallel if(_K0 > 100)
	{
#pragma omp for
		for (index_t _Row = 0; _Row < _K0; ++_Row) {
			_R0 = _R0 + pow(_Z, _Row)*_Data.rview(_Row);
		}
	}

	_Buff.rview(0) = _R0;
	for (size_t _Row = 1; _Row < _Height; ++_Row) {
		_Buff.rview(_Row) = _Data.rview(_Row) + _Buff.rview(_Row - 1)*_Z;
	}

	_Mycoeff.rview(_Height - 1) = _A * (_Buff.rview(_Height - 1) + _Buff.rview(_Height - 2)*_Z);
	for (index_t _Row = _Height - 2; _Row >= 0; --_Row) {
		_Mycoeff.rview(_Row) = _Z * (_Mycoeff.rview(_Row + 1) - _Buff.rview(_Row));
	}

	//Recursion over each row

	//pre-calculate the initial value for the first recursion
	matrix_type _C0(_Height, 1, zero<value_type>::value);
#pragma omp parallel if(_K0 > 100)
	{
#pragma omp for
		for (index_t _Col = 0; _Col < _K0; ++_Col) {
			_C0 = _C0 + _Mycoeff.cview(_Col) * pow(_Z, _Col);
		}
	}

	_Buff.cview(0) = _C0;
	for (std::size_t _Col = 1; _Col < _Width; ++_Col) {
		_Buff.cview(_Col) = _Mycoeff.cview(_Col) + _Buff.cview(_Col - 1)*_Z;
	}

	_Mycoeff.cview(_Width - 1) = _A * (_Buff.cview(_Width - 1) + _Buff.cview(_Width - 2)*_Z);
	for (index_t _Col = _Width - 2; _Col >= 0; --_Col) {
		_Mycoeff.cview(_Col) = _Z * (_Mycoeff.cview(_Col + 1) - _Buff.cview(_Col));
	}
}
template<typename _Ty>
void _Spline_interpolation<_Ty, _BIQNSPL>::_Coeff_impl() {
}
template<typename _Ty>
void _Spline_interpolation<_Ty, _BISPSPL>::_Coeff_impl() {
}

template class _Spline_interpolation<float, _BICBSPL>;
template class _Spline_interpolation<double, _BICBSPL>;
template class _Spline_interpolation<float, _BIQNSPL>;
template class _Spline_interpolation<double, _BIQNSPL>;
template class _Spline_interpolation<float, _BISPSPL>;
template class _Spline_interpolation<double, _BISPSPL>;
MATRICE_ALGS_END