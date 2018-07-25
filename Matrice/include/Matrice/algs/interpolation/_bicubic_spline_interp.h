#pragma once

#include "_base.h"

MATRICE_ALGS_BEGIN
template<typename _Ty, size_t _Options>
class BicubicSplineInterp : public InterpBase_<BicubicSplineInterp<_Ty, _Options>>
{
	using base_t = InterpBase_<BicubicSplineInterp<_Ty, _Options>>;
public:
	enum { Options = _Options };
	using value_t = _Ty;
	using matrix_t = types::Matrix<value_t>;

	MATRICE_GLOBAL_FINL BicubicSplineInterp(const matrix_t& _Data) {
		_Build_coeff(_Data);
	}

private:
	MATRICE_HOST_INL void _Build_coeff(const matrix_t& _Data)
	{
		const auto _Width = _Data.cols(), _Height = _Data.rows();
		m_coeff.create(_Height, _Width, zero<value_t>::value);

		//initialization
		const value_t _Z = -2. + sqrt(value_t(3));
		const value_t _A = _Z / (_Z*_Z - 1);
		const size_t _K0 = static_cast<size_t>(log(m_eps)/log(abs(_Z))) + 1;

		matrix_t _Buff(_Height, _Width);

		//recursion over each column
		//pre-calculate the initial value for the first recursion
		matrix_t _C0(1, _Width, zero<value_t>::value);
		for (size_t _Row = 0; _Row < _K0; ++_Row) {
			_C0 = _C0 + matrix_t(1, _Width, _Data[_Row])*pow(_Z, _Row);
		}

		transform(_C0.begin(), _C0.end(), _Buff.rbegin(0));
		for (size_t _Row = 1; _Row < _Height; ++_Row) {
			
		}


	}

	using base_t::m_coeff;
};
MATRICE_ALGS_END
