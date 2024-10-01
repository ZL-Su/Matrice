#include <functional>
#include <numeric>

#include "private/_plain_exp.hpp"
#include "private/_memory.h"
#include "../private/generic_fast_math.hpp"

#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
#include <mkl.h>
#endif

#define _EXPOP_EXPLICIT_INSTANTIATION(_Type, _Desc, _Name) \
template _Type* Exp::Op::_##_Desc##_##_Name<_Type>::operator()(\
int, _Type*, _Type*) const;

MATRICE_EXPR_BEGIN
template<typename _Ty>
_Ty* Exp::Op::_Mat_inv<_Ty>::operator()(int M, _Ty* Out, _Ty* In) const
{
	if (!In) return (Out);
	if (M == 2) { internal::_inv2x2m(In, Out); return (Out); };
	if (M == 3) { internal::_inv3x3m(In, Out); return (Out); };
	if (In != Out) privt::fill_mem(In, Out, M*M);

	if constexpr (type_bytes<_Ty>::value == 4) {
#ifdef MATRICE_MATH_KERNEL==MATRICE_USE_MKL
		LAPACKE_sgetri(/*LAPACK_ROW_MAJOR*/101, M, Out, M, nullptr);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports a kernel with preprocessor definition of MATRICE_MATH_KERNEL=MATRICE_USE_MKL.");
#endif
	}
	if constexpr (type_bytes<_Ty>::value == 8) {
#ifdef MATRICE_MATH_KERNEL==MATRICE_USE_MKL
		LAPACKE_dgetri(/*LAPACK_ROW_MAJOR*/101, M, Out, M, nullptr);
#else
		DGELOM_ERROR("Undefined math kernel, matrice supports a kernel with preprocessor definition of MATRICE_MATH_KERNEL=MATRICE_USE_MKL.");
#endif
	}

	return (Out);
}

_EXPOP_EXPLICIT_INSTANTIATION(float, Mat, inv)
_EXPOP_EXPLICIT_INSTANTIATION(double, Mat, inv)
MATRICE_EXPR_END