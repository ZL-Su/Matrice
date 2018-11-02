#pragma once
#include <tuple>
#include "../_type_traits.h"

DGE_MATRICE_BEGIN _DETAIL_BEGIN

template<typename _Ty, 
	typename = std::enable_if_t<is_floating_point_v<_Ty>>> 
	struct _Lapack_kernel_impl {};

_DETAIL_END DGE_MATRICE_END

#include "inl\_lnalge.inl"