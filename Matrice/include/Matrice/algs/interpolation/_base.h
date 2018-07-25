#pragma once

#include "../../core/matrix.h"
#include "../../core/vector.h"
#include "../../core/solver.h"

MATRICE_ALGS_BEGIN
enum {
	INTERP = 18,

	BSPLINE = 2307,
	BILINEAR = 2153,
	BICUBIC = 2156,
	BIQUINTIC = 2157,
	BISEPTIC = 2158
};
template<typename _Derived> class InterpBase_
{
	using derived_t = _Derived;
	using matrix_t = typename derived_t::matrix_t;
	using value_t = typename matrix_t::value_t;
public:
	InterpBase_() {}

protected:
	const value_t m_eps = 1.0e-7;
	matrix_t m_coeff;
};

MATRICE_ALGS_END
