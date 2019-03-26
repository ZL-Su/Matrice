/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once

#include "../../core/matrix.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty>
class _Error_analysis {
	using value_type = _Ty;
	using array_type = Matrix<_Ty>;
public:
	MATRICE_HOST_INL 
	static auto mean_bias(const array_type& _Vals, const value_type _Ref) {
		return (_Vals - _Ref).sum() / _Vals.size();
	}

};
_DETAIL_END
DGE_MATRICE_END