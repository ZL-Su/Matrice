/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include "../../util/_macros.h"

#ifndef MATRICE_MATH_KERNEL
#define MATRICE_MATH_KERNEL MATRICE_USE_NAT
#endif

#if MATRICE_MATH_KERNEL == MATRICE_USE_FKL
#include <fkl.h>
#endif
#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
#include <mkl.h>
#endif