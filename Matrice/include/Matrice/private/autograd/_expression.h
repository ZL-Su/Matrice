/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include "../_matrix_exp.hpp"

///<brief> construct computational expression system </brief>
///<method> exp::forward(): for expression evaluation </method>
///<method> exp::backward(): for derivatives computation based on chain rule </method>

DGE_MATRICE_BEGIN namespace exp {
_DETAIL_BEGIN

template<typename _Derived> class _Expression_base
{
	using _Myt = _Expression_base;
	using _Mydt = _Derived;
public:
	_Expression_base() {};
	~_Expression_base() {};

	/**
	 * \return value of an expresstion
	 */
	MATRICE_GLOBAL_INL auto forward() const {
		return (static_cast<const _Mydt*>(this)->operator());
	}
	/**
	 * \return gradient of an expresstion
	 */
	MATRICE_GLOBAL_INL auto backward() const {

	}
private:

};

_DETAIL_END
} DGE_MATRICE_END