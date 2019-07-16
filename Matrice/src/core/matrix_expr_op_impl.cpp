#include <functional>
#include <numeric>
#ifdef __use_mkl__
#include <mkl.h>
#else
#include <fkl.h>
#endif
#include "../../include/Matrice/private/_plain_exp.hpp"
#include "../../include/Matrice/private/_memory.h"
#include "../private/generic_fast_math.hpp"

#define _EXPOP_EXPLICIT_INSTANTIATION(_Type, _Desc, _Name) \
template _Type* Expr::Op::_##_Desc##_##_Name<_Type>::operator()(\
int, _Type*, _Type*) const;

MATRICE_NAMESPACE_EXPR_BEGIN
template<typename _Ty>
_Ty * Expr::Op::_Mat_inv<_Ty>::operator()(int M, _Ty * Out, _Ty * In) const
{
	if (!In) return (Out);
	if (M == 2) { privt::_inv2x2m(In, Out); return (Out); };
	if (M == 3) { privt::_inv3x3m(In, Out); return (Out); };
	if (In != Out) privt::fill_mem(In, Out, M*M);
	
	if constexpr (type_bytes<_Ty>::value == 4) 
#ifdef __use_mkl__
		LAPACKE_sgetri(LAPACK_ROW_MAJOR, M, (float*)Out, M, nullptr);
#else
		flapk::_sginv((fkl::sptr)Out, M);
#endif
	if constexpr (type_bytes<_Ty>::value == 8)
#ifdef __use_mkl__
		LAPACKE_dgetri(LAPACK_ROW_MAJOR, M, (double*)Out, M, nullptr);
#else
		flapk::_dginv((fkl::dptr)Out, M);
#endif
	return (Out);
}

_EXPOP_EXPLICIT_INSTANTIATION(float, Mat, inv)
_EXPOP_EXPLICIT_INSTANTIATION(double, Mat, inv)
MATRICE_NAMESPACE_EXPR_END