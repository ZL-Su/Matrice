/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include <exception>
#include "_macros.h"

DGE_MATRICE_BEGIN _DETAIL_BEGIN
/**
 * \record source code location
 * \Example: _Source_location _Loc{__func__, __FILE__, __LINE__}
 */
struct _Source_location {
	const char* _Func = nullptr;
	const char* _File = nullptr;
	long        _Line;
};
_DETAIL_END DGE_MATRICE_END