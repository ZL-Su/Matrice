/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#pragma once

#include "interpolation\_interpolation.h"

DGE_MATRICE_BEGIN
enum {
	bilinear = algs::INTERP | algs::BILINEAR,
	bcspline = algs::_BICBSPL,
	bqspline = algs::_BIQNSPL,
	bsspline = algs::_BISPSPL
}; //interpolation algorithms

/*******************************************************************
	              Unified Interpolation Interface
	    Copyright (c) : Zhilong (Dgelom) Su, since 31/Jul/2018
 ******************************************************************/
template<
	typename _Ty, 
	typename _Tag = _TAG bicspl_tag,
	MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
using interpolation = algs::_Interpolation_wrapper<_Ty, _Tag>;

DGE_MATRICE_END
