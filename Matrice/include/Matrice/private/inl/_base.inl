/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include "../_matrix_base.hpp"

DGE_MATRICE_BEGIN
_TYPES_BEGIN
template<typename _Derived, typename _Traits, typename _Type>
template<typename _Rhs>
auto& Base_<_Derived, _Traits, _Type>::inplace_sub(const _Rhs& _Right) {
	if constexpr (is_scalar_v<_Rhs>) {
		simd::Packet_<value_type, >
	}

}
_TYPES_END
DGE_MATRICE_END